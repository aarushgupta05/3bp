"""make_fig3.py - Fig. 3, cost ratio of fixed-step leapfrog against TSALF at
matched accuracy.

Two panels: matched on maximum relative energy error, and matched on
short-horizon position error. The plotted value is leapfrog cost divided by
TSALF cost, so above one means TSALF is the cheaper method.

Reads the table written by make_speed_table_ns5001.py and re-runs nothing.

PINNED POINTS. A matched-accuracy ratio only means something when both methods
were actually run at a setting that meets the target. When a method's coarsest
setting tested already beats the target, the table reports the cost of that
coarsest run, which is an upper bound on what the method would need, not a
measurement of what it needs. Those points are drawn as open markers here. On
IC6 in the energy panel every point is of this kind, because leapfrog at
dt = 0.08 is already at 0.00242 percent and TSALF at eta = 0.2 is already at
0.0168 percent, so the flat 1.53 ratio is just 1251 divided by 820, the ratio
of the two cheapest runs in the sweep. It is not evidence that TSALF is cheaper
on IC6, and it does not contradict Fig. 1: IC6 needs no refinement from either
method, so the matched-accuracy question does not bind there.

Drawn to the CJSJ template's figure specification: single column width, 8 point
Times New Roman, axis labels in words rather than symbols or abbreviations,
units in parentheses, no chart title, 300 dpi.

Usage:  python make_fig3.py [speed_table.txt] [out.png]
"""
import sys, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
    "font.size": 8, "axes.linewidth": 0.8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7, "axes.labelsize": 8,
})
COLUMN_IN = 3.4
fig, axs = plt.subplots(2, 1, figsize=(COLUMN_IN, COLUMN_IN * 1.48))
PANELS = [
    (axs[0], "energy", "Target maximum relative\nenergy error (percent)",
     "(a) Matched on energy error"),
    (axs[1], "rms", "Target short-horizon\nposition error",
     "(b) Matched on short-horizon position error"),
]
seen = {}
for ax, mode, ylab, title in PANELS:
    xs = [r for c in CFG for _, r, _ in data[c][mode]]
    if not xs: continue
    lo, hi = min(xs), max(xs)
    edge = hi * 2.6                                   # parking spot for unreachable
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
            ax.plot([edge], [t], marker="s", ms=4.6, mfc="none", mew=1.1,
                    color=COL[c], ls="none", zorder=4)
    ax.axvline(1.0, color="k", lw=0.9, ls=":", zorder=1)
    if any_un:
        ax.axvspan(hi * 1.6, edge * 1.5, color="0.92", zorder=0)
    ys = [t for c in CFG for t, _, _ in data[c][mode]] + [t for c in CFG for t in unre[c][mode]]
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlim(lo / 1.8, edge * 1.5)
    ax.set_ylim(max(ys) * 1.55, min(ys) / 1.55)
    ax.set_xlabel("Cost ratio, leapfrog divided by TSALF", labelpad=2)
    ax.set_ylabel(ylab, labelpad=2, fontsize=7.5)
    ax.set_title(title, fontsize=8, loc="left", pad=3)
    ax.grid(alpha=0.25, which="both", lw=0.4, zorder=0)
    ax.text(0.015, 0.035, "leapfrog cheaper", transform=ax.transAxes,
            fontsize=6, color="0.35")
    ax.text(0.985, 0.035, "TSALF cheaper", transform=ax.transAxes,
            fontsize=6, color="0.35", ha="right")

# one shared legend for both panels
h = [Line2D([], [], color=COL[c], lw=1.1, marker="o", ms=3.6,
            ls="--" if c == "IC1" else "-",
            label=c + (" (dashed: interpolated\nacross the ejection band)" if c == "IC1" else ""))
     for c in CFG]
if seen.get("pin"):
    h.append(Line2D([], [], color="0.4", marker="o", ms=4.4, mfc="white", mew=1.1,
                    ls="none", label="target already met at the\ncoarsest setting tested"))
if seen.get("un"):
    h.append(Line2D([], [], color="0.4", marker="s", ms=4.6, mfc="none", mew=1.1,
                    ls="none", label="leapfrog cannot reach this\ntarget while bounded"))
fig.legend(handles=h, fontsize=6, loc="lower center", ncol=3, framealpha=0.95,
           handlelength=1.8, labelspacing=0.35, columnspacing=1.0, borderpad=0.4,
           bbox_to_anchor=(0.5, 0.0))

fig.tight_layout(rect=[0, 0.155, 1, 1])
fig.subplots_adjust(hspace=0.55)
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
fig.savefig(OUT, dpi=300)
print("wrote", OUT)
for c in CFG:
    for mode in ("energy", "rms"):
        pin = sum(1 for _, _, p in data[c][mode] if p)
        print(f"  {c} {mode:<7}: {len(data[c][mode])} plotted, {pin} pinned, "
              f"{len(unre[c][mode])} unreachable by leapfrog")
