"""Calibrate a fixed Pivot burn-in independently of production measurements.

Every chain starts from the straight walk.  For each simulated length, several
fixed-seed pilot chains are followed long enough to contain a late stationary
region.  The first sustained plateau is located separately for the end-to-end
distance, acceptance indicator, and normalized free path.  Production starts
after two additional stability windows and never before twice the latest
onset.  The script aborts rather than silently accepting a cutoff too close to
the end of a pilot trace.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

from project_paths import (CONFIG, ROOT, THERMALIZATION, create_output_tree,
                           object_sha256)

SOURCE = ROOT / "src" / "pivot_saw.cpp"
OBSERVABLES = ("Re2", "accepted", "free_path_fraction")


def load_protocol() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def length_groups(protocol: dict) -> dict[int, str]:
    pivot = protocol["pivot"]
    groups = {int(n): "short" for n in pivot["short_lengths"]}
    groups.update({int(n): "intermediate" for n in pivot["intermediate_lengths"]})
    groups.update({int(n): "long" for n in pivot["long_attempts"]})
    if len(groups) != (len(pivot["short_lengths"])
                       + len(pivot["intermediate_lengths"])
                       + len(pivot["long_attempts"])):
        raise ValueError("Pivot length groups overlap")
    return dict(sorted(groups.items()))


def compile_pivot(directory: Path) -> Path:
    executable = directory / "pivot_saw.exe"
    subprocess.run(["g++", "-std=c++17", "-O3", str(SOURCE), "-o", str(executable)],
                   check=True)
    return executable


def pilot_seed(base: int, length: int, replica: int) -> int:
    return base + length * 10_007 + replica * 1_000_003


def run_pilot(executable: Path, length: int, group: str, replica: int,
              protocol: dict, overwrite: bool) -> Path:
    calibration = protocol["thermalization_calibration"]
    raw_dir = THERMALIZATION / "raw" / f"N_{length:06d}"
    raw_dir.mkdir(parents=True, exist_ok=True)
    final = raw_dir / f"replica_{replica:02d}.csv"
    if final.exists() and not overwrite:
        return final
    partial = final.with_suffix(".partial.csv")
    result = raw_dir / f"replica_{replica:02d}.result.partial.csv"
    partial.unlink(missing_ok=True)
    result.unlink(missing_ok=True)
    command = [
        str(executable), "--lengths", str(length),
        "--attempts", str(calibration["pilot_attempts"][group]),
        "--thermalization", "0",
        "--sample-every", str(calibration["record_every"][group]),
        "--seed", str(pilot_seed(calibration["seed_base"], length, replica)),
        "--out", str(result), "--series-out", str(partial),
        "--no-autocorrelation", "--quiet",
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    if not partial.exists():
        raise RuntimeError(f"Pilot trace was not created: {partial}")
    partial.replace(final)
    result.unlink(missing_ok=True)
    return final


def first_sustained_plateau(table: pd.DataFrame, observable: str,
                            calibration: dict) -> tuple[int, float, float]:
    """Return onset, late mean and uncertainty from independent replicas."""
    # Smooth within each replica only.  The ensemble columns remain independent
    # units for the uncertainty calculation.
    window = 5 if observable == "Re2" else 11
    # A trailing window avoids using future states to make an early time look
    # equilibrated.  The timestamp therefore marks the end of the local block.
    smooth = table.rolling(window, center=False, min_periods=window).mean()
    mean = smooth.mean(axis=1).to_numpy(float)
    times = smooth.index.to_numpy(int)

    late_start = int(math.floor((1.0 - calibration["late_fraction"]) * len(smooth)))
    replica_late_means = smooth.iloc[late_start:].mean(axis=0)
    late_mean = float(replica_late_means.mean())
    late_sem = float(replica_late_means.sem())
    late_curve_sd = float(smooth.iloc[late_start:].mean(axis=1).std(ddof=1))
    if not np.isfinite(late_mean) or late_mean == 0:
        raise RuntimeError(f"Invalid late reference for {observable}")

    if observable == "Re2":
        floor = calibration["relative_tolerance_Re2"] * abs(late_mean)
    elif observable == "accepted":
        floor = calibration["absolute_tolerance_acceptance"]
    else:
        floor = calibration["absolute_tolerance_free_path_fraction"]
    # Early chains can have a very large between-replica SEM.  Letting that SEM
    # enlarge the band would incorrectly classify the transient as equilibrium.
    # The allowed stationary scatter is instead estimated only from the late
    # ensemble curve, with a small observable-specific floor.
    tolerance = max(floor, calibration["sigma_tolerance"] * late_curve_sd)
    stable = np.abs(mean - late_mean) <= tolerance
    consecutive = int(calibration["consecutive_stable_bins"])
    for index in range(1, len(times) - consecutive):
        local_stable = bool(np.all(stable[index:index + consecutive]))
        tail_mean = float(np.nanmean(mean[index:]))
        tail_ok = abs(tail_mean - late_mean) <= tolerance
        if local_stable and tail_ok:
            return int(times[index]), late_mean, late_sem
    raise RuntimeError(
        f"No sustained {observable} plateau found; increase the pilot duration"
    )


def analyse_length(length: int, group: str, paths: list[Path],
                   protocol: dict) -> tuple[dict, pd.DataFrame]:
    frames = []
    for replica, path in enumerate(sorted(paths)):
        frame = pd.read_csv(path)
        required = {"attempt", *OBSERVABLES}
        if not required.issubset(frame.columns):
            raise RuntimeError(f"Incomplete pilot trace {path}")
        frame["replica"] = replica
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)
    calibration = protocol["thermalization_calibration"]
    onsets: dict[str, int] = {}
    late: dict[str, tuple[float, float]] = {}
    pivoted: dict[str, pd.DataFrame] = {}
    for observable in OBSERVABLES:
        table = raw.pivot(index="attempt", columns="replica", values=observable)
        onset, late_mean, late_sem = first_sustained_plateau(
            table, observable, calibration
        )
        onsets[observable] = onset
        late[observable] = (late_mean, late_sem)
        pivoted[observable] = table

    stride = int(calibration["record_every"][group])
    pilot_attempts = int(calibration["pilot_attempts"][group])
    latest_onset = max(onsets.values())
    stable_span = int(calibration["consecutive_stable_bins"]) * stride
    guarded_onset = latest_onset + calibration["safety_factor"] * stable_span
    burn = int(math.ceil(max(calibration["safety_factor"] * latest_onset,
                             guarded_onset) / stride) * stride)
    maximum = int(calibration["maximum_burn_fraction"] * pilot_attempts)
    if burn > maximum:
        raise RuntimeError(
            f"N={length}: proposed burn-in {burn} exceeds the validated limit "
            f"{maximum}; increase the pilot duration"
        )

    attempts = sorted(raw.attempt.unique())
    summary = pd.DataFrame({"N": length, "attempt": attempts})
    for observable, table in pivoted.items():
        summary[f"{observable}_mean"] = table.mean(axis=1).to_numpy()
        summary[f"{observable}_sem"] = table.sem(axis=1).to_numpy()
        summary[f"{observable}_late_mean"] = late[observable][0]
        summary[f"{observable}_late_sem"] = late[observable][1]
        summary[f"{observable}_onset"] = onsets[observable]
    summary["production_burn_in"] = burn
    summary["pilot_replicas"] = len(paths)

    row = {
        "N": length,
        "group": group,
        "pilot_attempts": pilot_attempts,
        "record_every": stride,
        "Re2_onset": onsets["Re2"],
        "acceptance_onset": onsets["accepted"],
        "free_path_fraction_onset": onsets["free_path_fraction"],
        "latest_onset": latest_onset,
        "stable_window_span": stable_span,
        "production_burn_in": burn,
        "safety_factor": calibration["safety_factor"],
    }
    return row, summary


def calibrate(executable: Path, jobs: int, overwrite: bool) -> pd.DataFrame:
    create_output_tree()
    protocol = load_protocol()
    groups = length_groups(protocol)
    calibration_definition = {
        "length_groups": groups,
        "thermalization_calibration": protocol["thermalization_calibration"],
    }
    definition_path = THERMALIZATION / "run_definition.json"
    definition = {
        "protocol_sha256": object_sha256(calibration_definition),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    }
    if definition_path.exists() and not overwrite:
        existing = json.loads(definition_path.read_text(encoding="utf-8"))
        if existing != definition:
            raise RuntimeError(
                "Calibration code or protocol changed. Use --overwrite to "
                "regenerate the pilot traces."
            )
    definition_path.write_text(json.dumps(definition, indent=2), encoding="utf-8")
    replicas = int(protocol["thermalization_calibration"]["replicas"])
    tasks = [(length, group, replica)
             for length, group in groups.items() for replica in range(replicas)]
    paths: dict[int, list[Path]] = {length: [] for length in groups}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        futures = {
            pool.submit(run_pilot, executable, length, group, replica,
                        protocol, overwrite): (length, replica)
            for length, group, replica in tasks
        }
        for completed, future in enumerate(as_completed(futures), 1):
            length, _ = futures[future]
            paths[length].append(future.result())
            if completed % 16 == 0 or completed == len(tasks):
                print(f"thermalization pilots: {completed}/{len(tasks)}", flush=True)

    rows, curves = [], []
    for length, group in groups.items():
        row, curve = analyse_length(length, group, paths[length], protocol)
        rows.append(row)
        curves.append(curve)
    cutoffs = pd.DataFrame(rows).sort_values("N")
    ensemble = pd.concat(curves, ignore_index=True).sort_values(["N", "attempt"])
    cutoffs.to_csv(THERMALIZATION / "cutoffs.csv", index=False)
    ensemble.to_csv(THERMALIZATION / "ensemble_relaxation.csv", index=False)
    manifest = {
        "description": "independent fixed-seed pilot calibration; production data are not used",
        "protocol_sha256": object_sha256(calibration_definition),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "observables": list(OBSERVABLES),
        "selection_rule": (
            "first sustained plateau for each observable; the latest onset plus "
            "two full stability windows (and at least twice the onset), rounded "
            "to the trace spacing"
        ),
        "replicas_per_length": replicas,
    }
    (THERMALIZATION / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"Wrote calibrated burn-in values for {len(cutoffs)} lengths")
    return cutoffs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    protocol = load_protocol()
    groups = length_groups(protocol)
    replicas = protocol["thermalization_calibration"]["replicas"]
    print(f"Thermalization calibration: {len(groups)} lengths x {replicas} replicas")
    if args.dry_run:
        return 0
    with tempfile.TemporaryDirectory(prefix="saw_calibration_build_") as build:
        calibrate(compile_pivot(Path(build)), args.jobs, args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
