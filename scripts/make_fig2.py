"""make_fig2.py - Fig. 2, cost against energy error for the binary-single systems.

Single panel in the same house style as Fig. 1: one colour per binary
separation, circles for fixed-step leapfrog, stars for TSALF, diamonds for
IAS15, crosses for runs in which a body was ejected. Cost is force evaluations
per simulated year.

Drawn to the CJSJ template's figure specification: single column width, 8 point
Times New Roman, axis labels in words rather than symbols or abbreviations,
units in parentheses, no chart title, 300 dpi.

IAS15's true energy error runs from about 3.9e-13% to 4.1e-12%, several
decades below everything else on the panel. Plotting it honestly would flatten
the rest of the figure, so each IAS15 point is drawn pinned to a floor just
below the lowest real value, with a downward arrow marking that the true value
is off the bottom. The script prints the true values so the caption can be
checked against them.

Reads the configuration sweep rather than re-running anything, so the figure cannot
drift from the data the paper describes in Table 2.

NOTE ON THE SAMPLING GRID. This figure sits on the configuration sweep's output grid,
n_samples = 5000, while every other number in the paper is on n_samples =
5001. The sweep was never re-run at 5001. The grid affects max |dE/E0| because
that is a maximum over output samples; it does not affect the bounded and
ejected verdicts, the force-evaluation counts, or the shape of the comparison,
which are what Fig. 2 is read for. See REPRODUCE.md.

Usage:  python make_fig2.py [config_sweep_newton.json] [out.png]
"""
import sys, os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "outputs", "config_sweep_newton.json")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "figures", "fig2_cost_energy_binary_single.png")

CFG = ["BS_0.05", "BS_0.10", "BS_0.20", "BS_0.40"]
COL = {"BS_0.05": "#1f77b4", "BS_0.10": "#ff7f0e", "BS_0.20": "#2ca02c", "BS_0.40": "#d62728"}
LEAPFROG = ["leapfrog_0.04", "leapfrog_0.02", "leapfrog_0.01", "leapfrog_0.005"]

# ---------------------------------------------------------------- read
D = json.load(open(SRC, encoding="utf-8"))
R = D["results"]
def get(cfg, meth):
    return next((r for r in R if r["config"] == cfg and r["method"] == meth
                 and "error" not in r and r.get("max_dE_pct") is not None), None)

data = {}
for c in CFG:
    lf = [(get(c, m)["fe_per_yr"], get(c, m)["max_dE_pct"], not get(c, m)["bounded"])
          for m in LEAPFROG if get(c, m)]
    ts = get(c, "tsalf"); ia = get(c, "ias15")
    data[c] = dict(lf=sorted(lf), ts=ts, ias=ia)

missing = [c for c in CFG if not data[c]["lf"] or not data[c]["ts"] or not data[c]["ias"]]
if missing:
    sys.exit(f"no data parsed for {missing} - check {SRC}")

# IAS15 is decades below the rest; pin it to a floor below the lowest real point
# and mark it with a downward arrow rather than rescaling the whole panel.
real = ([y for c in CFG for _, y, _ in data[c]["lf"]]
        + [data[c]["ts"]["max_dE_pct"] for c in CFG])
FLOOR = min(real) / 3.0
ARROW_TO = FLOOR / 2.5

# ---------------------------------------------------------------- draw
# CJSJ template figure specification
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "Nimbus Roman", "DejaVu Serif"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7, "axes.labelsize": 8,
})
COLUMN_IN = 3.4                      # CJSJ two-column text block, one column
fig, ax = plt.subplots(figsize=(COLUMN_IN, COLUMN_IN))

for c in CFG:
    col = COL[c]
    xs = [p[0] for p in data[c]["lf"]]
    ys = [p[1] for p in data[c]["lf"]]
    ax.plot(xs, ys, color=col, lw=1.1, zorder=2)
    EDGE = dict(markeredgecolor="white", markeredgewidth=0.6, linestyle="none")
    for x, y, ejected in data[c]["lf"]:
        if ejected:
            ax.plot(x, y, marker="X", color=col, markersize=5.5, zorder=4, **EDGE)
        else:
            ax.plot(x, y, marker="o", color=col, markersize=4.2, zorder=5, **EDGE)
    t = data[c]["ts"]
    ax.plot(t["fe_per_yr"], t["max_dE_pct"],
            marker="X" if not t["bounded"] else "*", color=col,
            markersize=5.5 if not t["bounded"] else 7.5, zorder=3, **EDGE)
    i = data[c]["ias"]
    ax.plot(i["fe_per_yr"], FLOOR, marker="D", color=col, markersize=4.6, zorder=3, **EDGE)
    ax.annotate("", xy=(i["fe_per_yr"], ARROW_TO), xytext=(i["fe_per_yr"], FLOOR / 1.25),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=1.1,
                                mutation_scale=9, shrinkA=0, shrinkB=0), zorder=3)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(ARROW_TO / 4.0, max(real) * 1.9)
ax.set_xlabel("Force evaluations per simulated year", labelpad=4)
ax.set_ylabel("Maximum relative energy error (percent)", labelpad=4)
ax.grid(True, which="major", color="0.85", lw=0.5)
ax.grid(True, which="minor", color="0.93", lw=0.35)
ax.set_axisbelow(True)

handles = [Line2D([0], [0], color=COL[c], lw=1.6,
                  label=c.split("_")[1] + " AU") for c in CFG]
handles += [
    Line2D([0], [0], color="0.25", marker="o", markersize=4.2, lw=1.1, label="Leapfrog"),
    Line2D([0], [0], color="0.25", marker="*", markersize=7.5, lw=0, label="TSALF"),
    Line2D([0], [0], color="0.25", marker="D", markersize=4.6, lw=0, label="IAS15"),
    Line2D([0], [0], color="0.25", marker="X", markersize=5.5, lw=0, label="Ejected"),
]
ax.legend(handles=handles, loc="upper right", fontsize=7, ncol=2,
          title="Binary separation", title_fontsize=7,
          framealpha=0.95, edgecolor="0.7", borderpad=0.4, labelspacing=0.3,
          handlelength=1.6, columnspacing=1.0)

fig.tight_layout()
fig.savefig(OUT, dpi=300)
print("wrote", OUT)
print(f"  n_samples grid of this data: {D.get('NS')}   force model: {D.get('force_model')}")
print(f"  IAS15 plotted at the floor {FLOOR:.3g}%; true values:")
for c in CFG:
    i = data[c]["ias"]
    print(f"    {c}: {i['max_dE_pct']:.3g}%  at {i['fe_per_yr']:.1f} fe/yr")
for c in CFG:
    ej = sum(1 for _, _, e in data[c]["lf"] if e) + (0 if data[c]["ts"]["bounded"] else 1)
    print(f"  {c}: {len(data[c]['lf'])} leapfrog, 1 tsalf, 1 ias15, {ej} ejected")
