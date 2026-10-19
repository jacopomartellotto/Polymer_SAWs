"""Aggregate the saved Simple-sampling/Pivot replicas and plot comparisons."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from comparison_data import load_grid
from estimate_nu import fit, power_law
from figure_style import (apply_style, heading,
                          error_points, save_figure, SIMPLE_COLOR, PIVOT_COLOR,
                          PAPER, INK)
from project_paths import ANALYSIS, CONFIG, FIGURES, PIVOT_SUMMARY, create_output_tree


def main():
    create_output_tree()
    grid = load_grid()
    summary = pd.read_csv(PIVOT_SUMMARY).set_index("N")
    render_comparison(grid, summary)


def render_comparison(grid, short):
    fit_protocol = json.loads(CONFIG.read_text(encoding="utf-8"))["analysis"]
    pivot = short.sort_index()
    assert pivot.index.is_unique
    assert len(pivot) == 43 and pivot.index.max() == 20000
    apply_style()
    blue, orange = SIMPLE_COLOR, PIVOT_COLOR
    x = grid.N.to_numpy()

    def save(fig, name):
        save_figure(fig, FIGURES, name)
        plt.close(fig)

    def draw_success(ax):
        error_points(ax, pivot.index, pivot.acceptance_mean,
                     pivot.acceptance_sem, marker="s", size=21,
                     color=orange, label="Pivot", zorder=2)
        error_points(ax, x, grid.survival, grid.survival_sem,
                     color=blue, label="Simple sampling", hollow=True, zorder=3)
        ax.set(xlabel="$N$", ylabel="Success probability", xscale="log",
               yscale="log")
        ax.grid(axis="y", alpha=.22)
        ax.legend(frameon=False, fontsize=9)

    # Full-range comparison, without a zoom inset.
    simple_re2 = grid.Re2_mean.to_numpy(float)
    simple_sem = grid.Re2_sem.to_numpy(float)

    # Draw the same curve_fit model and parameters used in Figure 5.
    # Curves are restricted to the fitted ranges.
    def scaling_curve(n, y, sem, points):
        n = np.asarray(n, dtype=float)
        result = fit(n, y, sem)
        curve_n = np.geomspace(n.min(), n.max(), points)
        return curve_n, power_law(curve_n, result["A"], result["nu"])

    pivot_fit = pivot.loc[pivot.index >= int(fit_protocol["pivot_N_min"])]
    pivot_curve_n, pivot_curve_re2 = scaling_curve(
        pivot_fit.index,
        pivot_fit.Re2_mean_mean,
        pivot_fit.Re2_mean_sem,
        600,
    )
    direct_mask = x >= int(fit_protocol["simple_N_min"])
    direct_curve_n, direct_curve_re2 = scaling_curve(
        x[direct_mask],
        simple_re2[direct_mask],
        simple_sem[direct_mask],
        200,
    )

    fig, ax = plt.subplots(figsize=(8.2, 5.8), layout="constrained")
    heading(fig, "End-to-end distance")
    error_points(ax, pivot.index, pivot.Re2_mean_mean,
                 pivot.Re2_mean_sem, marker="s", size=21,
                 color=orange, label="Pivot", zorder=2)
    error_points(ax, x, simple_re2, simple_sem, marker="o", size=32,
                 color=blue, label="Simple sampling", hollow=True, zorder=3)
    ax.plot(pivot_curve_n, pivot_curve_re2, color=INK, lw=1.05,
            label="Pivot fit", zorder=1)
    ax.plot(direct_curve_n, direct_curve_re2, color=blue, lw=.95,
            linestyle="--", label="Simple fit", zorder=1)
    ax.set(xlabel="$N$", ylabel=r"$\langle R_e^2\rangle$",
           xscale="log", yscale="log")
    ax.legend(loc="lower right", ncol=2, columnspacing=1.2,
              handlelength=2.2, fontsize=8.5)

    save(fig, "01_estensione")

    fig, ax = plt.subplots(figsize=(8.2, 5.8), layout="constrained")
    heading(fig, "Success probability")
    draw_success(ax)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    save(fig, "02_success_probability")

    # Three requested quantities; no extrapolated simple values at large N.
    specs = [("Re2", "Re2_mean", "Re2_sem", "Re2_mean"),
             ("Failure", "failure", "failure_sem", "failure"),
             ("Normalized free path", "mean_free_path_fraction", "mean_free_path_fraction_sem", "free_path_fraction_mean")]
    g = grid.set_index("N")
    numeric = pd.DataFrame(index=pivot.index)
    display = pd.DataFrame(index=pivot.index)
    for title, gmean, gsem, pbase in specs:
        for method, means, sems in (("Simple", g[gmean], g[gsem]),
                                    ("Pivot", pivot[pbase + "_mean"], pivot[pbase + "_sem"])):
            numeric[f"{title}_{method}_mean"] = means
            numeric[f"{title}_{method}_sem"] = sems
            display[f"{title} ({method})"] = [format_uncertainty(means.get(n, np.nan), sems.get(n, np.nan)) for n in display.index]
    numeric.index.name = display.index.name = "N"
    numeric.to_csv(ANALYSIS / "comparison_three_observables.csv")
    # Compact table image for the directly comparable lengths.
    fig, ax = plt.subplots(figsize=(16, 7.5), layout="constrained")
    ax.axis("off")
    rows = display.loc[x].reset_index()
    labels = ["N", "Simple\n⟨Re²⟩", "Pivot\n⟨Re²⟩", "Simple\nFailure", "Pivot\nFailure", "Simple\nNormalized free path", "Pivot\nNormalized free path"]
    table = ax.table(cellText=rows.values, colLabels=labels, loc="center", cellLoc="center",
                     colWidths=[.04,.15,.15,.15,.15,.18,.18])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (r,c), cell in table.get_celld().items():
        cell.set_height(.075 if r == 0 else .048)
        cell.set_edgecolor(PAPER)
        cell.set_facecolor("#E7E7E7" if r == 0 else ("#F4F4F4" if r % 2 else "white"))
    ax.set_title("Shared lengths · mean ± SEM", pad=12)
    save(fig, "04_comparison_table")
    print(f"Created English plots with {len(pivot)} Pivot lengths and "
          f"{len(grid)} Simple-sampling lengths.")
    print(display.loc[[4,30,20000]].to_string())


def format_uncertainty(value, error):
    if not np.isfinite(value) or not np.isfinite(error):
        return "—"
    decimals = max(0, 1 - int(np.floor(np.log10(error)))) if error > 0 else 6
    return f"{value:.{decimals}f} ± {error:.{decimals}f}"


if __name__ == "__main__":
    main()
