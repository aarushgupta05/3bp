# Reproducing the CJSJ paper

Every figure, table and quoted number in *"Structure over Resolution: An Operating Envelope for Stable Integration of the Chaotic Three-Body Problem"* (Gupta), and where it comes from.

Rewritten Sat 19 Sep 2026. It replaces the version in the full research tree, which is out of date in three ways that mattered: it states the wrong sampling convention, quotes superseded energy values, and names a core module under its old filename. Do not copy that one forward.

## Environment

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

No environment variable is needed. Pure Newtonian point masses, `eps = 0`, is the default force model in every script, because that is what every number in the paper uses. `SIMON_FORCE_MODEL=soft` switches to the softened model and writes to separate filenames, so a softened run can never overwrite a Newtonian one.

Two environments have produced the paper's numbers:

| | run of 12 Aug 2026 | independent check, 11 and 19 Sep 2026 |
|---|---|---|
| platform | Linux x86-64 | Linux aarch64 |
| python | 3.11.14 | 3.10.12 |
| numpy | 2.3.5 | 2.2.6 |
| rebound | 5.1.1 | 3.28.0 |

The IC1 frontier reproduces to every printed digit across both, on different CPU architectures and different major versions of REBOUND. The paper reports the x86-64 column; keep one column throughout rather than mixing. REBOUND 4 and above will not build on aarch64 because `whfast512` ships x86 AVX-512 assembly. IAS15 itself is unchanged between 3.x and 5.x.

**One exception, and it is worth knowing before a referee finds it.** Everything driven by this repository's own integrators reproduces exactly: every fixed-step leapfrog number, every TSALF number and every heuristic-adaptive number in `outputs/config_sweep_newton.json` is bit-identical across the two environments. The IAS15 rows are not. IAS15 chooses its own steps, so a different REBOUND build lands on a slightly different step sequence:

| | REBOUND 5.1.1 | REBOUND 3.28.0 |
|---|---|---|
| IC1 reference cost | 81,384 fe | 83,344 fe |
| IC4 reference cost | 80,000 fe | 80,968 fe |
| binary-single IAS15 energy error | 3.9e-13% to 4.1e-12% | 2.2e-13% to 2.3e-12% |

The cost figures move by at most 2.4%, which is invisible on the logarithmic axis of Fig. 2. The energy errors are all at the level of double-precision round-off, where the number is a property of the arithmetic rather than of the method, and both ranges are several decades below anything else on the panel, which is the claim the paper actually makes about IAS15. The range quoted in the Fig. 2 caption is the REBOUND 5.1.1 one.

## Sampling convention

**All energy errors are maxima over n_samples = 5001 on a T = 100 yr run, with the single documented exception of Fig. 2.**

This is not a free choice. The output grid is a zero-order hold: the state recorded at a sample time is the state at the first step boundary at or after it. When the grid spacing is not a whole multiple of the timestep, a fixed-step method picks up a phase error of up to one step, and adaptive methods do not, so the comparison is no longer fair.

Grid spacing is `T/(n_samples - 1)`:

| n_samples | spacing | commensurate with |
|---|---|---|
| 5001 | 0.02 | every dt in the sweep: 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125 |
| 2001 | 0.05 | only 0.01 and finer; not 0.08, 0.04 or 0.02 |
| 5000 | 0.020004... | none |

5001 is the only grid tested that divides every timestep in the sweep. An earlier version of this file recommended 2001 on the grounds that 0.05 was commensurate with every timestep tested. That is arithmetically false for dt = 0.08, 0.04 and 0.02, and the recommendation is withdrawn.

`max |dE/E0|` is a maximum over output samples, so for a bounded run with a sharp close encounter it depends on the sampling density:

| dt | ns=1001 | ns=2001 | ns=5001 |
|---|---|---|---|
| 0.005 | 0.06607% | 0.1497% | **2.524%** |
| 0.00125 | 0.002354% | 0.01136% | **0.1942%** |

**Bounded and ejected verdicts, max-r values and force-evaluation counts are identical at every sampling density.** The ejection band, which is the paper's actual result, does not move. Stability claims are stated in those terms for exactly this reason.

**The exception.** Fig. 2 sits on `outputs/config_sweep_newton.json`, which is the configuration sweep at n_samples = 5000. The sweep was never re-run at 5001. The grid affects the `max |dE/E0|` values plotted on the y axis; it does not affect the bounded and ejected verdicts, the force-evaluation counts on the x axis, or the ordering and separation of the four binary configurations, which are what Fig. 2 is read for. Closing this means re-running the sweep at 5001 and replacing the figure in the manuscript.

## Figures

| Paper | Script | Output |
|---|---|---|
| **Fig. 1** cost vs energy error, four analytical configurations | `scripts/make_fig1.py` | `figures/fig1_cost_energy_analytical.png` |
| **Fig. 2** cost vs energy error, binary-single | `scripts/make_fig2.py` | `figures/fig2_cost_energy_binary_single.png` |
| **Fig. 3** cost ratio at matched accuracy | `scripts/make_fig3.py` | `figures/fig3_cost_ratio_matched_accuracy.png` |
| **Fig. 4** Earth-Moon separation over 100 yr | `scripts/make_fig4.py` | `figures/fig4_earth_moon_separation.png` |

The committed PNGs are the images that appear in the manuscript. Re-running the scripts redraws them. On the environment above, Figs. 1 and 3 come back byte-identical, Fig. 4 comes back with identical content and a few bytes of difference from the matplotlib build, and Fig. 2 comes back with the same data and minor differences in axis limits and tick placement, because its committed image predates `make_fig2.py`. Regenerating Fig. 2 at n_samples = 5001 and replacing the manuscript image is the clean fix for both that and the grid exception above, and it is the one item in this repository that needs a change to the paper.

Figs. 1 and 3 read `outputs/speed_table_ns5001.txt` and re-run nothing, so neither can drift from the numbers the paper quotes. Fig. 2 reads `outputs/config_sweep_newton.json`. Fig. 4 reads `outputs/sun_earth_moon_run.json`.

Figs. 1, 2 and 3 in the manuscript were replaced on 19 Sep 2026 so that each has a committed generator. The versions in v1.7 and earlier had no generator anywhere: a search across every `.py` in the research tree that day found nothing matching their distinguishing features. Do not reintroduce the older images.

## Tables

| Paper | Script | Output |
|---|---|---|
| **Table 1** analytical ICs and lambda_fit | ICs: `src/simon_core.py::_ICS_RAW`, lambda_fit: `scripts/fit_lambda.py` | `outputs/lambda_fit.txt` |
| **Table 2** binary-single and Sun-Earth-Moon configurations | `scripts/config_sweep.py` | `outputs/config_sweep_newton.json` |
| **Table 3** diagnostic panel | definitions only, `src/simon_core.py::metric_panel` | none |
| **Table 4** IC1 fixed-step frontier | `scripts/ic1_frontier_5001.py` | `outputs/IC1_frontier_ns5001.txt` |

## Quoted numbers

| Claim | Source |
|---|---|
| IC1 ejection band: bounded at dt = 0.04 (3.24%), EJECT at 0.02 (355.6%, 156.89 AU) and 0.01 (2,944%), bounded again at 0.005 (2.52%) | `outputs/IC1_frontier_ns5001.txt` |
| TSALF bounded at every eta on IC1; ejects at eta = 0.2 on IC4 | `outputs/IC1_frontier_ns5001.txt`, `outputs/speed_table_ns5001.txt` |
| Sun-Earth-Moon: energy 0.007867%, Earth-Moon separation 1.98950 AU against a true 0.00257 AU, 774x | `outputs/metric_trap_table.txt` |
| Matched-accuracy cost ratios | `outputs/speed_table_ns5001.txt` |
| lambda_fit per configuration | `outputs/lambda_fit.txt` |

## Regenerating from scratch

Run in this order. `make_speed_table_ns5001.py` must come before Figs. 1 and 3, because they read the table it writes.

```bash
python scripts/ic1_frontier_5001.py          # Table 4                       (~10 s)
python scripts/fit_lambda.py                 # Table 1 lambda                (~3 s)
python scripts/make_speed_table_ns5001.py    # matched-accuracy table        (~25 s)
python scripts/make_fig1.py                  # Fig. 1
python scripts/make_fig3.py                  # Fig. 3
python scripts/make_fig2.py                  # Fig. 2
python scripts/make_fig4.py                  # Fig. 4 + metric-trap table
```

Each script takes its input and output paths from its own location, so any of them can be run from anywhere, and each prints the files it wrote. Total wall time on the aarch64 environment above is about 45 seconds.

`ic1_frontier_5001.py` re-runs four stored softened anchors at dt = 0.04 before it does anything else, and stops rather than printing a table if they do not reproduce:

| IC | max abs(dE/E0) |
|---|---|
| IC1 | 3.8564% |
| IC3 | 0.6718% |
| IC4 | 13.7722% |
| IC6 | 0.0006202% |

### Two inputs that are committed rather than regenerated

**`outputs/config_sweep_newton.json`**, the configuration sweep behind Table 2 and Fig. 2. `scripts/config_sweep.py` regenerates it, but it is not part of the sequence above, for two reasons. Its IAS15 rows are REBOUND-build dependent, as set out under Environment, and `make_speed_table_ns5001.py` reads one number out of it, the IAS15 reference cost printed in each configuration header, so regenerating the sweep changes two header lines in `speed_table_ns5001.txt`. No number the paper quotes is affected. If you do re-run it, run `make_speed_table_ns5001.py` again afterwards and expect those two header lines to move.

**`outputs/sun_earth_moon_run.json`**, the Sun-Earth-Moon run behind Fig. 4. `scripts/config_sweep.py` skips Sun-Earth-Moon by default, because regenerating it means refetching live state vectors from JPL Horizons, which would put the paper's headline number at the mercy of a network call and of whatever ephemeris Horizons is serving that day. Set `SIMON_SKIP_SEM=0` to refetch, and expect the numbers to move. The Sun-Earth-Moon configuration is integrated with the pure Newtonian point-mass force model in every run of the sweep, so these rows are the paper's force model whichever run they came from.

## A note on lambda_fit

The four lambda values in `src/simon_core.py::_ICS_RAW` are hardcoded literals carried over from earlier work, and no script regenerated them. `fit_lambda.py` is that missing script. It reports the unfloored least-squares slope of log(separation) against time for a shadow trajectory displaced by 1e-8, over the window 1e-7 < separation < 1e-1, with a 12-point robustness sweep, an exponential-versus-linear discrimination and a renormalised Benettin cross-check. Table 1 quotes its single-shadow column.

Two results from it are load-bearing:

- The stored IC6 value of **-0.0189 is not reproducible.** The regenerated value is **+0.0218**, and the entire sweep is positive. The legacy estimator floors its return at `max(slope, 0.02)` and so cannot emit a negative number at all.
- **IC6 separates linearly, not exponentially**, R-squared 0.998 against 0.924 for the exponential fit, the only configuration where linear wins. A divergence rate is not a meaningful descriptor there.

lambda_fit is used in one place only: it sets the averaging window for short-horizon RMS, `t <= min(2/max(lambda, 1e-3), T)`. It enters no force evaluation, step-size rule or energy diagnostic. That is what makes it defensible to call it a finite-time divergence rate rather than a converged Lyapunov exponent, and that sentence still needs to reach the paper's Methods section.

## What is in `src/`

`simon_core.py` holds the N-body harness, the initial conditions, the IAS15 reference and the diagnostic panel. `integrators_symplectic.py` holds TSALF. Run `python src/integrators_symplectic.py` for its own self-test: a Kepler e = 0.9 energy-boundedness check over 20 orbits, and a forward-then-reversed round trip, which is the property the implicit symmetric stepsize is supposed to buy.

The `simon_` prefix and the `SIMON_FORCE_MODEL` variable are historical. They are left alone deliberately: renaming touches every import, and it is not worth doing while a submission is open.
