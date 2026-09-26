"""make_fig1.py - Fig. 1, cost against energy error for the four analytical configurations.

Single panel, one colour per configuration, circles for fixed-step leapfrog and
stars for TSALF, crosses for runs in which a body was ejected. Cost is force
evaluations per simulated year, so fe / T with T = 100 yr.

Drawn to the CJSJ template's figure specification: single column width, 8 point
Times New Roman, axis labels in words rather than symbols or abbreviations,
units in parentheses, no chart title, 300 dpi. Markers carry a white edge so
that coincident points read as two markers rather than one; IC1 and IC3 at
dt = 0.08 differ by 1.86 against 1.80 percent and would otherwise hide each
other entirely.

Reads the matched-accuracy speed table rather than re-running anything, so the
figure cannot drift from the numbers the paper quotes.

Usage:  python make_fig1.py [speed_table_ns5001.txt] [out.png]
"""
import sys, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "outputs", "speed_table_ns5001.txt")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures", "fig1_cost_energy_analytical.png")

T = 100.0
CFG = ["IC1", "IC3", "IC4", "IC6"]
COL = {"IC1": "#1f77b4", "IC3": "#ff7f0e", "IC4": "#2ca02c", "IC6": "#d62728"}

# ---------------------------------------------------------------- parse
# config header:   IC1  (lambda~0.166)   ias15 reference cost ~81384 fe
# sweep lines:     leapfrog: dt=0.08: 1.86% / 1251fe | dt=0.02: 356% / 5001fe EJ | ...
POINT = re.compile(r"(?:dt|eta)=([\d.]+):\s*([\d.eE+-]+)%\s*/\s*(\d+)fe(\s*EJ)?")

data = {c: {"leapfrog": [], "tsalf": []} for c in CFG}
cfg = None
for line in open(SRC, encoding="utf-8"):
    m = re.match(r"^(IC\d)\s+\(lambda", line)
    if m:
        cfg = m.group(1) if m.group(1) in data else None
        continue
    if cfg is None:
        continue
    m = re.match(r"\s*(leapfrog|tsalf)\s*:", line)
    if not m:
        continue
    meth = m.group(1)
    for ctrl, err, fe, ej in POINT.findall(line):
        data[cfg][meth].append((int(fe) / T, float(err), bool(ej.strip())))

missing = [c for c in CFG if not data[c]["leapfrog"] or not data[c]["tsalf"]]
if missing:
    sys.exit(f"no data parsed for {missing} - check the speed table format")

# ---------------------------------------------------------------- draw
# CJSJ template figure specification
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.labelsize": 8,
})
COLUMN_IN = 3.4                      # CJSJ two-column text block, one column
fig, ax = plt.subplots(figsize=(COLUMN_IN, COLUMN_IN))

for c in CFG:
    for meth, mk in (("leapfrog", "o"), ("tsalf", "*")):
        pts = sorted(data[c][meth])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, color=COL[c], lw=1.1, zorder=2)
        for x, y, ejected in pts:
            if ejected:
                ax.plot(x, y, marker="X", color=COL[c], markersize=5.5,
                        markeredgecolor="white", markeredgewidth=0.6,
                        linestyle="none", zorder=4)
            else:
                ax.plot(x, y, marker=mk, color=COL[c],
                        markersize=6.5 if meth == "tsalf" else 4.2,
                        markeredgecolor="white", markeredgewidth=0.6,
                        linestyle="none", zorder=5 if meth == "leapfrog" else 3)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Force evaluations per simulated year", labelpad=4)
ax.set_ylabel("Maximum relative energy error (percent)", labelpad=4)
ax.grid(True, which="major", color="0.85", lw=0.5)
ax.grid(True, which="minor", color="0.93", lw=0.35)
ax.set_axisbelow(True)

handles = [Line2D([0], [0], color=COL[c], lw=1.6, label=c) for c in CFG]
handles += [
    Line2D([0], [0], color="0.25", marker="o", markersize=4.2, lw=1.1, label="Leapfrog"),
    Line2D([0], [0], color="0.25", marker="*", markersize=6.5, lw=1.1, label="TSALF"),
    Line2D([0], [0], color="0.25", marker="X", markersize=5.5, lw=0, label="Ejected"),
]
# key below the axes: at 8 point the box no longer fits inside without
# covering the IC1 ejection band, which is the point of the figure
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.155),
          fontsize=8, ncol=4, framealpha=0.95, edgecolor="0.7", borderpad=0.4,
          labelspacing=0.3, handlelength=1.6, columnspacing=1.0)

fig.tight_layout()
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print("wrote", OUT)
for c in CFG:
    n_lf, n_ts = len(data[c]["leapfrog"]), len(data[c]["tsalf"])
    ej = sum(1 for m in ("leapfrog", "tsalf") for p in data[c][m] if p[2])
    print(f"  {c}: {n_lf} leapfrog, {n_ts} tsalf, {ej} ejected")
