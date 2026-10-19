"""Execute thermalization calibration, production, analysis, and plotting."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys

from project_paths import ROOT, create_output_tree


def run(command: list[str]) -> None:
    print("\n> " + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    create_output_tree()
    if shutil.which("g++") is None:
        raise RuntimeError("g++ was not found on PATH; a C++17 compiler is required")

    python = sys.executable
    common = ["--jobs", str(max(1, args.jobs))]
    if args.overwrite:
        common.append("--overwrite")
    if args.dry_run:
        run([python, "-B", str(ROOT / "scripts" / "calibrate_thermalization.py"),
             *common, "--dry-run"])
        run([python, "-B", str(ROOT / "scripts" / "run_production.py"),
             "all", *common, "--dry-run"])
        return 0

    run([python, "-B", str(ROOT / "scripts" / "calibrate_thermalization.py"),
         *common])
    run([python, "-B", str(ROOT / "scripts" / "run_production.py"),
         "all", *common])

    run([python, "-B", str(ROOT / "scripts" / "render_outputs.py")])
    print("\nPipeline complete. Final plots are in outputs/figures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
