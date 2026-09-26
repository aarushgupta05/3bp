"""make_fig3.py - Fig. 3, cost ratio of fixed-step leapfrog against TSALF at
matched accuracy.

Two panels: matched on maximum relative energy error, and matched on
short-horizon position error. The plotted value is leapfrog cost divided by
TSALF cost, so above one means TSALF is the cheaper method.

Reads the table written by make_speed_table_ns5001.py and re-runs nothing.

UNITS. Panel (a) is matched on maximum relative energy error, which is a
percentage. Panel (b) is matched on short-horizon position RMS, which is a
distance in the dimensionless units of Table I (G = 1), not a percentage.

PINNED POINTS. A matched-accuracy ratio only means something when both methods
were actually run at a setting that meets the target. When a method's coarsest
setting tested already beats the target, the table reports the cost of that
coarsest run, which is an upper bound on what the method would need, not a
measurement of what it needs. Those points are drawn as open circles here. On
IC6 in the energy panel every point is of this kind, because leapfrog at
dt = 0.08 is already at 0.00242 percent, below every target in the panel, so
its cost is pinned at 1251 force evaluations throughout. At the four coarsest
targets TSALF is likewise pinned, at eta = 0.2 and 820 force evaluations, so
the 1.53 ratio there is just 1251 divided by 820, the ratio of the two cheapest
runs in the sweep. It is not evidence that TSALF is cheaper on IC6, and it does
not contradict Fig. 1: IC6 needs no refinement from either method, so the
matched-accuracy question does not bind there.

UNREACHABLE TARGETS. Some targets fixed-step leapfrog never meets while every
body stays bound, however fine the timestep, so the cost ratio is not a finite
number and there is nowhere on the axis to put it. Those points sit in a
narrow panel past an axis break, under the tick "never". The break marks say
the axis is cut there, so the squares are not read as very large ratios.

Drawn to the CJSJ template's figure specification: single column width, 8 point
Times New Roman throughout including the legend and every in-plot note, axis
labels in words rather than symbols or abbreviations, units in parentheses, no
chart title, 300 dpi.

Usage:  python make_fig3.py [speed_table.txt] [out.png]
"""
import sys, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "outputs", "speed_table_ns5001.txt")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures", "fig3_cost_ratio_matched_accuracy.png")
CFG = ["IC1", "IC3", "IC4", "IC6"]
COL = {"IC1": "#1f77b4", "IC3": "#ff7f0e", "IC4": "#2ca02c", "IC6": "#d62728"}

# ---------------------------------------------------------------- parse
TEXT = open(SRC, encoding="utf-8").read()

def cheapest_bounded(block, key):
    """Cost of the cheapest run of `key` that stayed bounded, in force evals."""
    line = re.search(rf"{key}\s*: (.*)", block).group(1)
    best = None
    for seg in line.split("|"):
        m = re.search(r"(?:dt|eta)=([\d.]+):\s*([\d.eE+-]+)%\s*/\s*(\d+)fe(\s*EJ)?", seg)
        if m and not m.group(4):
            fe = int(m.group(3))
            best = fe if best is None or fe < best else best
    return best

FLOOR = {}
for c in CFG:
    blk = re.search(rf"^{c}\s+\(lambda.*?(?=^={{10,}}|\Z)", TEXT, re.S | re.M).group(0)
    FLOOR[c] = (cheapest_bounded(blk, "leapfrog"), cheapest_bounded(blk, "tsalf"))

cfg = mode = None
data = {c: {"energy": [], "rms": []} for c in CFG}   # (target, ratio, pinned)
unre = {c: {"energy": [], "rms": []} for c in CFG}
for line in TEXT.splitlines():
    m = re.match(r"^(IC\d)\s+\(lambda", line)
    if m: cfg = m.group(1); continue
    if "target max|dE|" in line: mode = "energy"; continue
    if "target rms_short" in line: mode = "rms"; continue
    if line.startswith("====="): mode = None; continue
    if not (cfg and mode): continue
    m = re.match(r"\s+([\d.]+)%?\s+(\d+)\s+(\d+)\s+([\d.]+)x", line)
    if m:
        lf_fe, ts_fe = int(m.group(2)), int(m.group(3))
        pinned = (lf_fe == FLOOR[cfg][0]) or (ts_fe == FLOOR[cfg][1])
        data[cfg][mode].append((float(m.group(1)), float(m.group(4)), pinned))
        continue
    m = re.match(r"\s+([\d.]+)%?\s+--\s+\d+\s+LF n/a", line)
    if m:
        unre[cfg][mode].append(float(m.group(1)))

# ---------------------------------------------------------------- draw
PT = 8
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
    "font.size": PT, "axes.linewidth": 0.8,
    "xtick.labelsize": PT, "ytick.labelsize": PT,
    "legend.fontsize": PT, "axes.labelsize": PT, "axes.titlesize": PT,
    "mathtext.fontset": "custom", "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
})
COLUMN_IN = 3.4
fig = plt.figure(figsize=(COLUMN_IN, 4.55))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 0.13], wspace=0.06,
              hspace=0.52, left=0.215, right=0.975, top=0.960, bottom=0.215)

PANELS = [
    (0, "energy", "Target maximum relative\nenergy error (percent)",
     "(a) Matched on energy error",
     [1.0, 0.3, 0.1, 0.03, 0.01], ["1.0", "0.3", "0.1", "0.03", "0.01"],
     [0.3, 1, 3, 10], ["0.3", "1", "3", "10"]),
    (1, "rms", "Target short-horizon position\nerror (dimensionless units)",
     "(b) Matched on position error",
     [0.05, 0.02, 0.01, 0.005, 0.002], ["0.05", "0.02", "0.01", "0.005", "0.002"],
     [0.3, 0.5, 0.7, 1.0], ["0.3", "0.5", "0.7", "1.0"]),
]
seen = {}
for row, mode, ylab, title, yticks, ylabels, xticks, xlabels in PANELS:
    ax = fig.add_subplot(gs[row, 0])
    bx = fig.add_subplot(gs[row, 1], sharey=ax)
    xs = [r for c in CFG for _, r, _ in data[c][mode]]
    lo, hi = min(xs), max(xs)
    any_un = any(unre[c][mode] for c in CFG)

    for c in CFG:
        pts = sorted(data[c][mode])
        if pts:
            y = [p[0] for p in pts]; x = [p[1] for p in pts]
            ax.plot(x, y, "--" if c == "IC1" else "-", lw=1.1, color=COL[c], zorder=3)
            for t, r, pinned in pts:
                if pinned:
                    seen["pin"] = True
                    ax.plot([r], [t], marker="o", ms=4.4, mfc="white", mew=1.1,
                            color=COL[c], ls="none", zorder=5)
                else:
                    ax.plot([r], [t], marker="o", ms=3.6, color=COL[c],
                            mec="white", mew=0.5, ls="none", zorder=4)
        for t in unre[c][mode]:
            seen["un"] = True
            bx.plot([0.5], [t], marker="s", ms=4.6, mfc="none", mew=1.1,
                    color=COL[c], ls="none", zorder=4, clip_on=False)

    ax.axvline(1.0, color="k", lw=0.9, ls=":", zorder=1)
    ys = [t for c in CFG for t, _, _ in data[c][mode]] + [t for c in CFG for t in unre[c][mode]]
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlim(lo / 1.9, hi * 1.9)
    ax.set_ylim(max(ys) * 1.7, min(ys) / 1.7)
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels)
    ax.set_yticks([], minor=True)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels)
    ax.set_xticks([], minor=True)
    ax.set_ylabel(ylab, labelpad=2)
    ax.set_title(title, loc="left", pad=3)
    ax.grid(alpha=0.25, which="both", lw=0.4, zorder=0)
    ax.text(0.015, 0.030, "leapfrog cheaper", transform=ax.transAxes, color="0.35")
    ax.text(0.985, 0.030, "TSALF cheaper", transform=ax.transAxes,
            color="0.35", ha="right")
    # x label centred under the pair of axes
    ax.set_xlabel("Cost ratio, leapfrog divided by TSALF", labelpad=2, x=0.565)

    # the narrow off-scale panel, past an axis break
    bx.set_xlim(0, 1)
    bx.set_xticks([0.5]); bx.set_xticklabels(["never"])
    bx.tick_params(axis="y", which="both", left=False, labelleft=False)
    bx.grid(axis="y", alpha=0.25, which="both", lw=0.4, zorder=0)
    ax.spines["right"].set_visible(False)
    bx.spines["left"].set_visible(False)
    if not any_un:
        bx.set_visible(False)
        ax.spines["right"].set_visible(True)
        ax.set_xlabel("Cost ratio, leapfrog divided by TSALF", labelpad=2, x=0.5)
        continue
    # break marks on the facing spines
    kw = dict(marker=[(-1, -0.6), (1, 0.6)], ms=5, ls="none", mec="k", mew=0.9,
              color="k", clip_on=False)
    ax.plot([1, 1], [0, 1], transform=ax.transAxes, **kw)
    bx.plot([0, 0], [0, 1], transform=bx.transAxes, **kw)

# one shared legend for both panels
h = [Line2D([], [], color=COL[c], lw=1.1, marker="o", ms=3.6,
            ls="--" if c == "IC1" else "-", label=c) for c in CFG]
if seen.get("pin"):
    h.append(Line2D([], [], color="0.4", marker="o", ms=4.4, mfc="white", mew=1.1,
                    ls="none", label="target already met at the coarsest setting"))
fig.legend(handles=h[:4], loc="lower center", ncol=4, frameon=False,
           handlelength=1.9, columnspacing=1.0, handletextpad=0.5,
           bbox_to_anchor=(0.5, 0.058))
fig.legend(handles=h[4:], loc="lower center", ncol=1, frameon=False,
           handlelength=1.9, handletextpad=0.5, bbox_to_anchor=(0.5, 0.002))

os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
fig.savefig(OUT, dpi=300)
print("wrote", OUT)
for c in CFG:
    for mode in ("energy", "rms"):
        pin = sum(1 for _, _, p in data[c][mode] if p)
        print(f"  {c} {mode:<7}: {len(data[c][mode])} plotted, {pin} pinned, "
              f"{len(unre[c][mode])} unreachable by leapfrog")
