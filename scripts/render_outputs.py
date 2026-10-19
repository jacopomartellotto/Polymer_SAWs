"""Recalculate the analysis and render the final thesis figures."""
import subprocess
import sys

from project_paths import ANALYSIS, FIGURES, ROOT, THERMALIZATION

FIGURE_NAMES = (
    "01_estensione", "02_success_probability", "04_comparison_table",
    "05_nu_scaling", "07_thermalization", "09_nu_vs_lmin",
    "backup_pivot_corrections",
)


def main() -> None:
    scripts = (
        "compare_simple_sampling_pivot.py",
        "estimate_nu.py",
        "plot_nu_scaling.py",
        "plot_thermalization.py",
        "plot_nu_vs_lmin.py",
        "plot_corrected_fit.py",
    )
    for name in scripts:
        print(f"Running {name}", flush=True)
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / name)],
                       cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL if name == "estimate_nu.py" else None)

    required = [FIGURES / f"{name}.png" for name in FIGURE_NAMES]
    required += [
        ANALYSIS / "comparison_three_observables.csv",
        ANALYSIS / "nu_estimates.json",
        ANALYSIS / "05_nu_scaling_metadata.json",
        ANALYSIS / "nu_vs_lmin.csv",
        ROOT / "slides" / "sampling_fit_results.tex",
        THERMALIZATION / "cutoffs.csv",
        THERMALIZATION / "ensemble_relaxation.csv",
    ]
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing pipeline outputs:\n" + "\n".join(map(str, missing)))
    print(f"Validated {len(required)} outputs")


if __name__ == "__main__":
    main()
