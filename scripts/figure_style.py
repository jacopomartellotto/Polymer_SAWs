"""Standard Matplotlib figures with supplied SEMs and explicit zoom regions."""
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import LogLocator, NullFormatter

PAPER = PANEL = INSET = "#FFFFFF"
INK = "#222222"
SIMPLE_COLOR, PIVOT_COLOR = "#1f77b4", "#ff7f0e"
MUTED = "#222222"


def apply_style():
    plt.style.use("default")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11, "mathtext.fontset": "dejavusans",
        "figure.facecolor": PAPER, "axes.facecolor": PANEL,
        "savefig.facecolor": PAPER, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.edgecolor": INK, "axes.linewidth": .8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.labelsize": 12, "axes.labelpad": 5,
        "axes.titlesize": 10, "axes.titleweight": "normal", "axes.titlepad": 7,
        "axes.grid": True, "axes.axisbelow": True, "grid.alpha": .18,
        "grid.color": "#B0B0B0", "grid.linewidth": .7,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.top": False, "ytick.right": False,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "xtick.major.width": .55, "ytick.major.width": .55,
        "xtick.minor.width": .4, "ytick.minor.width": .4,
        "legend.frameon": False, "legend.fontsize": 10,
        "legend.handlelength": 1.5, "legend.labelspacing": .3,
        "legend.handletextpad": .5, "legend.borderaxespad": .7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


def error_points(ax, x, y, sem, *, color, marker="o", label=None,
                 size=24, hollow=False, zorder=3):
    """Draw every measurement with its precomputed, unchanged SEM."""
    xlabel, ylabel = ax.get_xlabel(), ax.get_ylabel()
    ax.errorbar(x, y, yerr=sem, fmt="none", color=color,
                elinewidth=.75, capsize=1.8, capthick=.65, zorder=zorder-1)
    ax.scatter(x, y, marker=marker, s=size*(.7 if hollow else .6),
               facecolors="none" if hollow else color,
               edgecolors=color, linewidths=.8 if hollow else .4,
               label=label, zorder=zorder)
    # Series names are internal columns, not publication axis labels.
    ax.set(xlabel=xlabel, ylabel=ylabel)


def heading(fig, title, footer="43 Pivot lengths · 16 Simple-sampling lengths · error bars: ±1 SEM"):
    # Paper figures use external captions, not decorative titles or footers.
    fig.set_label(title)
    fig.set_size_inches((12, 5.2) if len(fig.axes) > 1 else (7.2, 5.2))


def inset_style(ax):
    ax.set_facecolor(INSET)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#888888")
        spine.set_linewidth(.7)
    ax.tick_params(labelsize=8, pad=2)
    ax.title.set_fontsize(9)
    ax.title.set_color(INK)


def connect_zoom(ax, inset, direction="above", external=False):
    """Outline the visible zoom bounds and link facing corners in data coordinates.

    Clip the source rectangle to the parent view, including a positive lower
    bound when a linear inset includes zero on a logarithmic parent axis.
    Artists do not participate in autoscaling or constrained layout.
    """
    x0, x1 = inset.get_xlim()
    y0, y1 = inset.get_ylim()
    px0, px1 = ax.get_xlim()
    py0, py1 = ax.get_ylim()
    x0, x1 = max(x0, px0), min(x1, px1)
    y0, y1 = max(y0, py0), min(y1, py1)
    if x1 <= x0 or y1 <= y0:
        raise ValueError("Zoom window does not intersect the parent view")
    style = dict(edgecolor="#9AA4AE", linewidth=.45, zorder=1.1)
    rectangle = Rectangle((x0, y0), x1-x0, y1-y0, fill=False,
                          transform=ax.transData, **style)
    rectangle.set_in_layout(False)
    ax.add_artist(rectangle)
    if direction == "above":
        corners = [((x0, y1), (0, 0)), ((x1, y1), (1, 0))]
    elif direction == "right":
        corners = [((x1, y0), (0, 0)), ((x1, y1), (0, 1))]
    else:
        raise ValueError(f"Unknown zoom direction: {direction}")
    for source, target in corners:
        line = ConnectionPatch(source, target, coordsA=ax.transData,
                               coordsB=inset.transAxes, clip_on=False,
                               **style)
        line.set_in_layout(False)
        (ax.figure if external else ax).add_artist(line)


def panel(ax, letter, title):
    ax.annotate(f"({letter.lower()})", xy=(0, 1), xycoords="axes fraction",
                xytext=(0, 5), textcoords="offset points",
                ha="left", va="bottom", fontsize=11, color=INK,
                annotation_clip=False, zorder=20)


def save_figure(fig, out, name):
    def clean(ax, child=False):
        ax.grid(False, which="minor")
        ax.grid(True, which="major", alpha=.14 if child else .18)
        ax.tick_params(which="both", labelsize=8 if child else 10,
                       direction="in" if child else "out", top=False, right=False)
        if child:
            ax.tick_params(which="major", length=2.3, width=.45)
            ax.tick_params(which="minor", length=1.3, width=.35)
        for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
            if scale == "log":
                axis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
                axis.set_minor_formatter(NullFormatter())
        legend = ax.get_legend()
        if legend:
            for text in legend.get_texts():
                text.set_fontsize(10)
        for child_ax in ax.child_axes:
            clean(child_ax, True)
    for ax in fig.axes:
        clean(ax)
    fig.savefig(out / f"{name}.png",
                dpi=600 if name != "04_comparison_table" else 220)
