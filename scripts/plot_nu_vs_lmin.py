"""Plot the Pivot exponent estimate against the minimum fitted length."""
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.optimize import curve_fit

from estimate_nu import fit
from figure_style import apply_style, save_figure, PIVOT_COLOR, INK
from project_paths import ANALYSIS, FIGURES

TABLE = ANALYSIS / "nu_vs_lmin.csv"


def calculate() -> pd.DataFrame:
    source = pd.read_csv(ANALYSIS / "comparison_three_observables.csv")
    data = source.dropna(subset=["Re2_Pivot_mean", "Re2_Pivot_sem"]).copy()
    # The pure-power-law stability study starts at N=200 and retains at least
    # four measured lengths (two residual degrees of freedom) per fit.
    data = data[data.N >= 200].sort_values("N")
    rows = []
    # Every simulated length is tested as a cutoff while at least four measured
    # lengths remain.  N_min is selected; nu is fitted from all N >= N_min.
    for lmin in data.N:
        subset = data[data.N >= lmin]
        if len(subset) < 4:
            continue
        result = fit(subset.N, subset.Re2_Pivot_mean,
                     subset.Re2_Pivot_sem)
        rows.append({
            # Keep the legacy CSV field name for compatibility with the saved
            # analysis, while using N_min in figures and explanatory text.
            "L_min": int(lmin), "L_max": int(subset.N.max()),
            "nu": result["nu"], "nu_ci68_low": result["ci68"][0],
            "nu_ci68_high": result["ci68"][1],
            "chi2": result["chi2"], "dof": result["dof"],
            "nu_sem": result["nu_sem"], "chi2_reduced": result["chi2_reduced"],
            "p_value": result["p_value"], "method": result["method"],
            "absolute_sigma": result["absolute_sigma"],
            "points": result["points"], "model": "A*N^(2*nu)",
        })
    table = pd.DataFrame(rows)
    table.to_csv(TABLE, index=False)
    return table


def select_constant_fit(table: pd.DataFrame) -> dict:
    """Try the last 9, 8, ... points, stopping at reduced chi-square <=1.5.

    User-requested descriptive diagonal fit: shared-data correlations between
    cutoff estimates are ignored. Its nominal SEM and chi-square are not a
    calibrated plateau test and do not include range-selection uncertainty.
    At least three points are required (two residual degrees of freedom).
    """
    ordered = table.sort_values("L_min")
    trials = []
    selected = None
    for count in range(min(9, len(ordered)), 2, -1):
        subset = ordered.tail(count)
        x = subset.L_min.to_numpy(float)
        y = subset.nu.to_numpy(float)
        sem = subset.nu_sem.to_numpy(float)
        popt, pcov = curve_fit(
            lambda n, c: np.full_like(n, c), x, y,
            sigma=sem, absolute_sigma=True, p0=[.75],
            jac=lambda n, c: np.ones((len(n), 1)),
        )
        constant = float(popt[0])
        chi2 = float(np.sum(((y-constant)/sem)**2))
        result = dict(points=count, N_min=int(x.min()), N_max=int(x.max()),
                      constant=constant, constant_sem=float(np.sqrt(pcov[0, 0])),
                      chi2=chi2, dof=count-1, chi2_reduced=chi2/(count-1))
        trials.append(result)
        if result["chi2_reduced"] <= 1.5:
            selected = result
            break
    record = dict(selected=selected, trials=trials,
                  rule="Last 9, then 8, ... cutoff estimates; first reduced chi-square <=1.5; minimum 3 points",
                  method="curve_fit; sigma=nu_sem; absolute_sigma=True",
                  caveat="Diagonal errors: correlations ignored. Nominal uncertainty conditional on selected range; descriptive selection, not a calibrated plateau-onset test.")
    (ANALYSIS / "nu_constant_region.json").write_text(json.dumps(record, indent=2))
    if selected is None:
        raise ValueError("No tested constant interval meets the 1.5 threshold")
    print(pd.DataFrame(trials)[["points", "N_min", "chi2_reduced"]].to_string(index=False))
    return selected


def plot(table: pd.DataFrame) -> None:
    """Display the cutoff estimates and nominal 1/2-sigma bands of the constant."""
    apply_style()
    constant_fit = select_constant_fit(table)
    # Display starts at 500; retain the complete cutoff scan in the CSV.
    # No estimate is invented at 500: the first measured cutoff shown is 600.
    table = table[table.L_min >= 500].copy()
    x = table.L_min.to_numpy(float)
    y = table.nu.to_numpy(float)
    low = table.nu_ci68_low.to_numpy(float)
    high = table.nu_ci68_high.to_numpy(float)
    selected_cutoff = constant_fit["N_min"]
    selected = table[table.L_min == selected_cutoff]
    if len(selected) != 1:
        raise ValueError("Configured Pivot cutoff is missing from the stability scan")
    fig, ax = plt.subplots(figsize=(7.6, 5.35), layout="constrained")
    # Mark the operational constant-fit selection, not a proven plateau onset.
    accent = "#355C73"
    ax.axvspan(selected_cutoff, x.max()*1.06, color=accent, alpha=.04,
               linewidth=0, zorder=0)
    ax.axvline(selected_cutoff, color=accent, lw=1.0, ls=(0, (3, 3)), zorder=2)
    constant = constant_fit["constant"]
    constant_sem = constant_fit["constant_sem"]
    # These are conditional Gaussian intervals for the fitted constant;
    # they are not prediction bands for individual cutoff estimates.
    ax.fill_between([selected_cutoff, x.max()], constant - 2*constant_sem,
                    constant + 2*constant_sem, color="#89BDD3", alpha=.32,
                    linewidth=0, zorder=1)
    ax.fill_between([selected_cutoff, x.max()], constant - constant_sem,
                    constant + constant_sem, color="#3F88A8", alpha=.36,
                    linewidth=0, zorder=1.1)
    ax.errorbar(x, y, yerr=np.vstack((y - low, high - y)),
                fmt="s", linestyle="none", markersize=4.2,
                color=PIVOT_COLOR, markeredgecolor="white",
                markeredgewidth=.5, elinewidth=.8, capsize=2,
                capthick=.7, zorder=3)
    ax.scatter([selected_cutoff], selected.nu, s=52, marker="s",
               facecolors="none", edgecolors=accent, linewidths=1.1, zorder=4)
    theory_line = ax.axhline(.75, color=INK, lw=.95, ls=(0, (4, 2.5)),
                             label=r"Theory: $\nu=3/4$")
    constant_line = ax.hlines(constant, selected_cutoff, x.max(),
                              color=accent, lw=1.5, label="Constant fit")
    ax.legend([theory_line, constant_line,
               Patch(facecolor="#3F88A8", alpha=.36),
               Patch(facecolor="#89BDD3", alpha=.32)],
              [r"Theory: $\nu=3/4$", "Constant fit",
               r"$\pm1\sigma$ on constant", r"$\pm2\sigma$ on constant"],
              loc="upper left", fontsize=8.6)
    # The chi-square criterion applies to entire candidate constant-fit ranges,
    # rather than to an individual cutoff estimate or its underlying power fit.
    ax.text(.38, .70,
            "Excluded from constant fit\n"
            + r"($\chi^2_{\mathrm{red}}$ too large)",
            transform=ax.transAxes, ha="center", va="center", fontsize=9.0,
            color="#57636D", linespacing=1.45)
    decimals = max(0, -int(np.floor(np.log10(constant_fit["constant_sem"]))))
    summary = (
        rf"$\nu={constant_fit['constant']:.{decimals}f}\pm{constant_fit['constant_sem']:.{decimals}f}$"
        + "\n" + rf"$\chi^2/\mathrm{{ndof}}={constant_fit['chi2_reduced']:.2f}$"
        + "\n" + rf"$N_{{\min}}\geq {selected_cutoff}$"
    )
    ax.text(.975, .055, summary,
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9.8,
            color=INK, linespacing=1.7,
            bbox=dict(boxstyle="round,pad=.8,rounding_size=.5",
                      facecolor="white", edgecolor="#94AFC0", linewidth=.8))
    margin = max(3e-5, .10 * (high.max() - low.min()))
    ax.set(ylim=(min(low.min(), .75) - margin,
                 max(high.max(), .75) + margin),
           xlabel=r"$N_{\min}$", ylabel=r"$\nu$")
    ax.set_xscale("log")
    ax.set_xlim(500, x.max()*1.06)
    ticks = sorted({500, 1000, 2000, selected_cutoff, 10000})
    ax.set_xticks(ticks, [str(v) for v in ticks])
    fig.set_label("Pivot exponent stability: cutoff fits with one-sigma covariance errors")
    save_figure(fig, FIGURES, "09_nu_vs_lmin")
    plt.close(fig)


def main() -> None:
    table = calculate()
    plot(table)
    print(f"Pivot: {len(table)} cutoff fits from L_min={table.L_min.min()} "
          f"to {table.L_min.max()}; first fit uses {table.points.max()} lengths.")


if __name__ == "__main__":
    main()
