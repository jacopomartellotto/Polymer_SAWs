"""Plot end-to-end-distance scaling and the fitted exponent."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from estimate_nu import fit, power_law, format_nu_estimate
from comparison_data import load_grid
from figure_style import (apply_style, heading, error_points, save_figure,
                          SIMPLE_COLOR, PIVOT_COLOR, INK)
from project_paths import ANALYSIS, CONFIG, FIGURES


def main():
    protocol = json.loads(CONFIG.read_text(encoding="utf-8"))["analysis"]
    pivot_nmin = int(protocol["pivot_N_min"])
    simple_nmin = int(protocol["simple_N_min"])
    pivot_nmin_tex = f"{pivot_nmin:,}".replace(",", r"\,")
    table = pd.read_csv(ANALYSIS / "comparison_three_observables.csv")
    simple = table.dropna(subset=["Re2_Simple_mean"])
    fit_data = table[table.N >= pivot_nmin]
    n = fit_data.N.to_numpy(float)
    y = fit_data.Re2_Pivot_mean.to_numpy()
    sem = fit_data.Re2_Pivot_sem.to_numpy()
    result = fit(n, y, sem)
    nu = result["nu"]
    curve_n = np.geomspace(n.min(), n.max(), 600)
    curve_re2 = power_law(curve_n, result["A"], result["nu"])
    simple_data = simple[simple.N >= simple_nmin]
    simple_result = fit(simple_data.N, simple_data.Re2_Simple_mean,
                        simple_data.Re2_Simple_sem)
    simple_n = simple_data.N.to_numpy(float)
    simple_curve_n = np.geomspace(simple_n.min(), simple_n.max(), 200)
    simple_curve_re2 = power_law(simple_curve_n, simple_result["A"],
                                simple_result["nu"])
    blue, orange, dark = SIMPLE_COLOR, PIVOT_COLOR, INK
    apply_style()
    fig, ax = plt.subplots(figsize=(8.2, 5.8), layout="constrained")
    heading(fig, "The scaling exponent")
    fig.set_size_inches(8.2, 5.8)
    grid_counts = load_grid().set_index("N").survivors
    reliable = simple.N.map(grid_counts).to_numpy() >= 100

    def points(ax, power, subset=table, include_simple=True, labels=True):
        scale = subset.N.to_numpy(float)**power
        error_points(ax, subset.N, subset.Re2_Pivot_mean / scale,
                     subset.Re2_Pivot_sem / scale, marker="s", color=orange, size=21,
                     label="Pivot" if labels else None, zorder=3)
        if include_simple:
            scale = simple.N.to_numpy(float)**power
            for mask, marker, label in ((reliable, "o", "Simple sampling"),
                                        (~reliable, "x", "Simple sampling (<100 samples)")):
                if not mask.any():
                    continue
                error_points(ax, simple.N[mask], simple.Re2_Simple_mean[mask] / scale[mask],
                             simple.Re2_Simple_sem[mask] / scale[mask], marker=marker,
                             color=blue, hollow=True,
                             label=label if labels else None, zorder=4)

    points(ax, 0)
    ax.plot(curve_n, curve_re2,
            color=dark, lw=1, label="Pivot fit", zorder=2)
    ax.plot(simple_curve_n, simple_curve_re2,
            color=blue, lw=.9, linestyle="--", label="Simple fit", zorder=2)
    ax.set(xscale="log", yscale="log", xlabel="$N$",
           ylabel=r"$\langle R_e^2\rangle$")
    ax.grid(axis="y", alpha=.22)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    uncertainty = result["nu_sem"]
    simple_uncertainty = simple_result["nu_sem"]
    annotation = (
        "Pivot: " + r"$\langle R_e^2\rangle=A N^{2\nu}$" + "\n"
        + r"$\nu=" + format_nu_estimate(nu, uncertainty) + "$; "
        + rf"${pivot_nmin_tex}\leq N\leq20\,000$" + "\n\n"
        + "Simple sampling: " + r"$\langle R_e^2\rangle=A N^{2\nu}$" + "\n"
        + r"$\nu=" + format_nu_estimate(simple_result['nu'], simple_uncertainty) + "$; "
        + rf"${simple_nmin}\leq N\leq30$" + "\n\n"
        + r"Theory: $\nu=3/4$")
    ax.text(.035, .965, annotation, transform=ax.transAxes,
            va="top", fontsize=9, linespacing=1.45,
            bbox={"facecolor":"white", "edgecolor":"none", "alpha":.9, "pad":3})
    save_figure(fig, FIGURES, "05_nu_scaling")
    plt.close(fig)
    result["plotted_pivot_points"] = len(table)
    result["simple_fit"] = simple_result
    result["notes"] = (
        f"Nu is fitted with the pure law A*N^(2*nu), using Pivot N>={pivot_nmin} "
        f"and Simple sampling N>={simple_nmin}. Curves are shown only in their fit ranges. "
        "The Simple-sampling exponent is an effective short-chain estimate. "
        "Error bars are one SEM. Nu errors are sqrt(pcov[1,1]) from curve_fit with absolute_sigma=True, conditional on the model and fit range. "
        "Distances use lattice spacing 1."
    )
    (ANALYSIS / "05_nu_scaling_metadata.json").write_text(json.dumps(result, indent=2))
    print(f"Plotted {len(table)} Pivot and {len(simple)} Simple-sampling lengths; free nu={nu:.9f}.")


if __name__ == "__main__":
    main()
