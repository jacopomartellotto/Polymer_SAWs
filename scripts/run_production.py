"""Run the calibrated thesis simulations and aggregate independent replicas."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

import pandas as pd

from project_paths import (CONFIG, PIVOT, PIVOT_RUNS, PIVOT_SUMMARY, ROOT,
                           SIMPLE, THERMALIZATION, create_output_tree,
                           object_sha256)

SOURCE = ROOT / "src" / "pivot_saw.cpp"
SUMMARY_COLUMNS = ["Re2_mean", "acceptance", "failure", "free_path_mean",
                   "free_path_fraction_mean"]


def load_protocol() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def compile_pivot(directory: Path) -> Path:
    executable = directory / "pivot_saw.exe"
    subprocess.run(["g++", "-std=c++17", "-O3", str(SOURCE), "-o", str(executable)],
                   check=True)
    return executable


def pivot_jobs(protocol: dict, cutoffs: pd.DataFrame) -> list[dict]:
    pivot = protocol["pivot"]
    cutoff_by_n = cutoffs.set_index("N").production_burn_in.to_dict()
    jobs = []
    groups = (
        (pivot["short_lengths"], pivot["short_attempts"],
         pivot["replicas_short_and_intermediate"], "short"),
        (pivot["intermediate_lengths"], pivot["intermediate_attempts"],
         pivot["replicas_short_and_intermediate"], "intermediate"),
        ([int(n) for n in pivot["long_attempts"]], None,
         pivot["replicas_long"], "long"),
    )
    for lengths, common_attempts, replicas, group in groups:
        for length in lengths:
            length = int(length)
            if length not in cutoff_by_n:
                raise RuntimeError(f"Missing calibrated burn-in for N={length}")
            attempts = (int(pivot["long_attempts"][str(length)])
                        if group == "long" else int(common_attempts))
            for replica in range(1, int(replicas) + 1):
                jobs.append({
                    "N": length,
                    "group": group,
                    "replica": replica,
                    "attempts": attempts,
                    "thermalization": int(cutoff_by_n[length]),
                })
    return jobs


def result_path(job: dict) -> Path:
    return (PIVOT_RUNS / f"N_{job['N']:06d}"
            / f"replica_{job['replica']:04d}.csv")


def validate_existing(path: Path, job: dict) -> None:
    frame = pd.read_csv(path)
    if len(frame) != 1:
        raise RuntimeError(f"Expected one result row in {path}")
    row = frame.iloc[0]
    expected = {"N": job["N"], "attempts": job["attempts"],
                "thermalization": job["thermalization"]}
    for key, value in expected.items():
        if int(row[key]) != int(value):
            raise RuntimeError(
                f"Existing {path} has {key}={row[key]}, expected {value}; "
                "rerun with --overwrite"
            )


def run_one(executable: Path, job: dict, protocol: dict, overwrite: bool) -> Path:
    final = result_path(job)
    if final.exists() and not overwrite:
        validate_existing(final, job)
        return final
    final.parent.mkdir(parents=True, exist_ok=True)
    partial = final.with_suffix(".partial.csv")
    partial.unlink(missing_ok=True)
    seed = (int(protocol["pivot"]["production_seed_base"])
            + job["replica"] * 1_000_003 + job["N"] * 1_000_000_007)
    command = [
        str(executable), "--lengths", str(job["N"]),
        "--attempts", str(job["attempts"]),
        "--thermalization", str(job["thermalization"]),
        "--sample-every", "1", "--seed", str(seed),
        "--out", str(partial), "--no-autocorrelation", "--quiet",
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    partial.replace(final)
    return final


def aggregate(paths: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted(paths)]
    raw = pd.concat(frames, ignore_index=True)
    rows = []
    for length, frame in raw.groupby("N", sort=True):
        row = {
            "N": int(length),
            "replicas": len(frame),
            "attempts_per_replica": int(frame.attempts.iloc[0]),
            "thermalization": int(frame.thermalization.iloc[0]),
        }
        for metric in SUMMARY_COLUMNS:
            values = frame[metric].astype(float)
            row[f"{metric}_mean"] = values.mean()
            row[f"{metric}_sd"] = values.std(ddof=1)
            row[f"{metric}_sem"] = values.sem(ddof=1)
        rows.append(row)
    summary = pd.DataFrame(rows).sort_values("N")
    PIVOT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(PIVOT_SUMMARY, index=False)
    return summary


def check_reuse_manifest(protocol: dict, overwrite: bool) -> None:
    manifest_path = PIVOT / "production_manifest.json"
    if not manifest_path.exists() or overwrite:
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        "protocol_sha256": object_sha256(protocol["pivot"]),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    }
    if any(manifest.get(key) != value for key, value in expected.items()):
        raise RuntimeError(
            "Code or protocol changed after the existing Pivot runs were made. "
            "Use --overwrite to regenerate them."
        )


def run_pivot(executable: Path, protocol: dict, jobs_count: int,
              overwrite: bool) -> None:
    cutoff_path = THERMALIZATION / "cutoffs.csv"
    if not cutoff_path.exists():
        raise RuntimeError("Run calibrate_thermalization.py before production")
    cutoffs = pd.read_csv(cutoff_path)
    jobs = pivot_jobs(protocol, cutoffs)
    check_reuse_manifest(protocol, overwrite)
    definition_path = PIVOT / "run_definition.json"
    definition = {
        "protocol_sha256": object_sha256(protocol["pivot"]),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "cutoffs_sha256": hashlib.sha256(cutoff_path.read_bytes()).hexdigest(),
    }
    if definition_path.exists() and not overwrite:
        existing = json.loads(definition_path.read_text(encoding="utf-8"))
        if existing != definition:
            raise RuntimeError(
                "Production code, protocol, or calibrated cutoffs changed. "
                "Use --overwrite to regenerate the Pivot runs."
            )
    PIVOT.mkdir(parents=True, exist_ok=True)
    definition_path.write_text(json.dumps(definition, indent=2), encoding="utf-8")
    paths = []
    with ThreadPoolExecutor(max_workers=max(1, jobs_count)) as pool:
        futures = {
            pool.submit(run_one, executable, job, protocol, overwrite): job
            for job in jobs
        }
        for completed, future in enumerate(as_completed(futures), 1):
            paths.append(future.result())
            if completed % 32 == 0 or completed == len(jobs):
                print(f"Pivot production: {completed}/{len(jobs)} jobs", flush=True)
    summary = aggregate(paths)
    manifest = {
        "description": "production starts after independently calibrated fixed burn-in",
        "protocol_sha256": object_sha256(protocol["pivot"]),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "cutoffs_sha256": hashlib.sha256(cutoff_path.read_bytes()).hexdigest(),
        "seed_formula": "production_seed_base + replica*1000003 + N*1000000007",
        "lengths": summary.N.astype(int).tolist(),
        "raw_results": len(paths),
    }
    (PIVOT / "production_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"Wrote {PIVOT_SUMMARY}")


def run_simple(protocol: dict, overwrite: bool) -> None:
    simple = protocol["simple_sampling"]
    manifest_path = SIMPLE / "production_manifest.json"
    definition = {
        "protocol_sha256": object_sha256(simple),
        "source_sha256": hashlib.sha256(
            (ROOT / "grid" / "run_simple_sampling.py").read_bytes()
        ).hexdigest(),
    }
    if manifest_path.exists() and not overwrite:
        saved = json.loads(manifest_path.read_text(encoding="utf-8"))
        if any(saved.get(key) != value for key, value in definition.items()):
            raise RuntimeError(
                "Simple-sampling code or protocol changed. Use --overwrite to regenerate it."
            )
        if (SIMPLE / "refined" / "grid_refined.csv").exists():
            print("Simple-sampling data already match the protocol; reusing them")
            return
    lengths = ",".join(map(str, simple["lengths"]))
    command = [
        sys.executable, str(ROOT / "grid" / "run_simple_sampling.py"),
        "--lengths", lengths,
        "--walks", str(simple["base_walks_per_length"]),
        "--seed", str(simple["base_seed"]),
        "--out", str(SIMPLE / "base.csv"),
        "--refinement-out", str(SIMPLE / "refined"),
        "--protocol", str(CONFIG),
        "--thesis",
    ]
    if overwrite:
        command.append("--force-refinement")
    started = time.perf_counter()
    subprocess.run(command, cwd=ROOT, check=True)
    definition.update(
        proposal_scheme=simple["proposal_scheme"],
        attempts_total=(len(simple["lengths"]) * int(simple["base_walks_per_length"])
                        + sum(map(int, simple["refinement_walks"].values()))),
        seconds=time.perf_counter() - started,
        completed_at_utc=datetime.now(timezone.utc).isoformat(),
    )
    manifest_path.write_text(json.dumps(definition, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", choices=["simple-sampling", "pivot", "all"])
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    create_output_tree()
    protocol = load_protocol()
    if args.selection in {"pivot", "all"}:
        cutoff_path = THERMALIZATION / "cutoffs.csv"
        if cutoff_path.exists():
            print(f"Pivot production jobs: {len(pivot_jobs(protocol, pd.read_csv(cutoff_path)))}")
        else:
            print("Pivot production requires thermalization calibration")
    if args.selection in {"simple-sampling", "all"}:
        print("Simple sampling: base run plus fixed high-statistics refinement")
    if args.dry_run:
        return 0
    if args.selection in {"simple-sampling", "all"}:
        run_simple(protocol, args.overwrite)
    if args.selection in {"pivot", "all"}:
        with tempfile.TemporaryDirectory(prefix="saw_production_build_") as build:
            run_pivot(compile_pivot(Path(build)), protocol, args.jobs, args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
