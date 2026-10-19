"""Plot the independently calibrated Pivot thermalization protocol."""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from figure_style import apply_style, panel, save_figure, PIVOT_COLOR, SIMPLE_COLOR, INK
from project_paths import FIGURES, THERMALIZATION

# Plot thermalization figures of some representive lenghts
REPRESENTATIVE_LENGTHS = (30, 150, 1_000, 20_000)

#Plot function
def main() -> None:
    data = pd.read_csv(THERMALIZATION / "ensemble_relaxation.csv")
    cutoffs = pd.read_csv(THERMALIZATION / "cutoffs.csv").set_index("N")
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), layout="constrained")
    for letter, ax, length in zip("ABCD", axes.flat, REPRESENTATIVE_LENGTHS):
        local = data[data.N == length].sort_values("attempt")
        row = cutoffs.loc[length]
        onset = int(row.latest_onset)
        burn = int(row.production_burn_in)
        reference = float(local.Re2_late_mean.iloc[0])
        reference_sem = float(local.Re2_late_sem.iloc[0])
        t = local.attempt.to_numpy()
        mean = local.Re2_mean.to_numpy() / reference
        sem = local.Re2_sem.to_numpy() / reference
        visible = t <= min(1.2 * burn, t.max())

        ax.axvspan(0, onset, color="#F3E6E3", alpha=.75, linewidth=0,
                   label="Detected transient" if length == 30 else None)
        ax.axhspan(max(1 - reference_sem / reference, 1e-6),
                   1 + reference_sem / reference, color=PIVOT_COLOR,
                   alpha=.14, linewidth=0,
                   label="Late plateau" if length == 30 else None)
        ax.fill_between(t[visible], np.maximum(mean[visible] - sem[visible], 1e-6),
                        mean[visible] + sem[visible], color=PIVOT_COLOR,
                        alpha=.18, linewidth=0)
        ax.plot(t[visible], mean[visible], color=PIVOT_COLOR, lw=1.05,
                label=r"Ensemble $\langle R_e^2(t)\rangle$" if length == 30 else None)
        ax.axhline(1, color=INK, lw=.75, ls="--")
        ax.axvline(onset, color=SIMPLE_COLOR, lw=1.1,
                   label="Latest plateau onset" if length == 30 else None)
        ax.axvline(burn, color=INK, lw=1.15, ls=":",
                   label="Production burn-in" if length == 30 else None)
        ax.annotate(f"{burn:,}".replace(",", "\u2009"),
                    xy=(burn, 0), xycoords=("data", "axes fraction"),
                    xytext=(-4, 6), textcoords="offset points", rotation=90,
                    ha="right", va="bottom", fontsize=8, color=INK)
        ax.set(xlim=(0, min(1.2 * burn, t.max())), yscale="log",
               xlabel="Attempted pivots, $t$",
               ylabel=r"$\langle R_e^2(t)\rangle/\langle R_e^2\rangle_{\rm late}$")
        ax.set_title(f"$N={length:,}$".replace(",", r"\,"))
        panel(ax, letter, "")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=5,
               handlelength=1.7, columnspacing=1.0, fontsize=8.5)
    fig.set_label("Independent Pivot burn-in calibration")
    save_figure(fig, FIGURES, "07_thermalization")
    plt.close(fig)


if __name__ == "__main__":
    main()
