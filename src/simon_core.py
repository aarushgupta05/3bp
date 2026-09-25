"""
simon_core.py - the N-body harness behind every number in the paper.

One integrator core, general over N bodies, plus the initial conditions, the
IAS15 reference and the diagnostic panel. Everything in scripts/ imports it.

What is here:
  * Initial conditions IC1, IC3, IC4 and IC6, centre-of-mass centred. IC2 and
    IC5 are defined but the paper does not use them.
  * A pairwise force model with two settings, selected by mode:
      newton : unsoftened Newtonian for every pair. This is the model the
               paper uses throughout, with eps = 0 in the energy diagnostic.
      soft   : Plummer-softened for pairs closer than SOFT_THRESH, Newtonian
               beyond it. Retained only so the four softened anchors checked
               by scripts/ic1_frontier_5001.py still reproduce.
  * Fixed-step leapfrog, with optional heuristic adaptive sub-stepping.
  * A REBOUND IAS15 reference trajectory.
  * The diagnostic panel of Table 3: max |dE/E0|, short-horizon position RMS
    against IAS15, force-evaluation cost, and a bounded or ejected verdict.

Sampling. The output grid is index-derived, t_cur = (step + 1) * dt, so there
is no accumulation drift. Any sample slot left unwritten is filled with the
final state, and two assertions catch an unwritten or all-zero slot rather
than letting it reach a table. All energy errors in the paper are maxima over
n_samples = 5001 on a T = 100 yr run.

Conventions. IC1 to IC6 use G = 1.0. The Sun-Earth-Moon system uses
G = 4 * pi^2, with distances in AU and time in years.
"""
import math
import numpy as np

# ---- physical / model constants ----
EPS          = 3e-4          # Plummer softening length, used by mode='soft' only
SOFT_THRESH  = 500.0 * EPS   # 0.15 AU: softening applies only to closer pairs
R_SOFT_MIN   = 5e-4          # below this softened separation, fall back to Newtonian
ADAPT_THRESH = 0.05          # default adaptive sub-step trigger
MAX_SUBSTEPS = 16            # default adaptive cap
G_TOY  = 1.0
G_REAL = 4.0 * np.pi**2

EJECT_AU = 10.0              # bounded := max dist from COM < this

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
_ICS_RAW = {
    "IC1": dict(m=[1.0, 0.01, 0.005], x0=[[0,0,0],[1,0,0],[0,1.2,0]],
                v0=[[0,0,0],[0,1,0],[-0.9,0,0]], exp_lambda=0.1655, stable=True),
    "IC2": dict(m=[1.0, 0.5, 0.25], x0=[[0,0,0],[1,0,0],[-0.5,0.8,0]],
                v0=[[0,0,0],[0,0.6,0],[-0.4,-0.3,0]], exp_lambda=0.0430, stable=False),
    "IC3": dict(m=[1.0, 0.01, 0.005], x0=[[0,0,0],[0.5,0,0],[0,2.5,0]],
                v0=[[0,0,0],[0,1.3,0],[-0.4,0,0]], exp_lambda=0.1207, stable=True),
    "IC4": dict(m=[1.0, 0.01, 0.005], x0=[[0,0,0],[1,0,0],[0,5.0,0]],
                v0=[[0,0,0],[0,1.0,0],[-0.12,0,0]], exp_lambda=0.0990, stable=True),
    "IC5": dict(m=[1.0, 0.01, 0.005], x0=[[0,0,0],[0.3,0,0],[0,2.0,0]],
                v0=[[0,0,0],[0,2.2,0],[-0.3,0,0]], exp_lambda=0.0648, stable=False),
    "IC6": dict(m=[1.0, 0.01, 0.005], x0=[[0,0,0],[1.5,0,0],[-3.0,0,0]],
                v0=[[0,0,0],[0,0.816,0],[0,-0.577,0]], exp_lambda=-0.0189, stable=True),
}
# The exp_lambda values above are hardcoded literals carried over from earlier
# work; no script regenerated them. scripts/fit_lambda.py is that script, and
# Table 1 of the paper quotes its unfloored least-squares values, not these.
# The stored IC6 value of -0.0189 is not reproducible: the regenerated value is
# +0.0218, and IC6 separates linearly rather than exponentially. lambda enters
# one place only, the averaging window for short-horizon RMS in metric_panel.

def get_ic(name):
    """Return (m, x0, v0) COM-centred, float64."""
    d = _ICS_RAW[name]
    m  = np.array(d["m"], dtype=np.float64)
    x0 = np.array(d["x0"], dtype=np.float64)
    v0 = np.array(d["v0"], dtype=np.float64)
    M = m.sum()
    x0 = x0 - (m[:, None] * x0).sum(0) / M
    v0 = v0 - (m[:, None] * v0).sum(0) / M
    return m, x0, v0

def ic_meta(name):
    return dict(_ICS_RAW[name])

# ---------------------------------------------------------------------------
# Force model. mode in {'newton','soft'} controls close-pair handling:
#   newton : pure unsoftened Newtonian for all pairs (what the paper uses)
#   soft   : Plummer-softened for pairs closer than SOFT_THRESH, Newtonian far
# Returns acc (N,3). Force-eval accounting is the caller's job, 1 per call.
# ---------------------------------------------------------------------------
def make_pairs(N):
    ii, jj = [], []
    for i in range(N):
        for j in range(i+1, N):
            ii.append(i); jj.append(j)
    return np.array(ii), np.array(jj)

def compute_acc(x, ii, jj, Gmimj, inv_mi, inv_mj, eps2, mode):
    rij = x[jj] - x[ii]
    r2  = np.einsum('ij,ij->i', rij, rij)
    r   = np.sqrt(r2 + 1e-30)
    invr3 = 1.0 / (r2 * r + 1e-30)
    F = Gmimj * invr3                       # Newtonian default (all pairs)
    if mode != 'newton':
        close = r < SOFT_THRESH
        if np.any(close):
            r2c = r2[close]
            r_soft_c = np.sqrt(r2c + eps2)
            denom = (r2c + eps2) ** 1.5 + 1e-30
            F_soft = Gmimj[close] / denom
            fb = (r_soft_c < R_SOFT_MIN)
            F[close] = np.where(fb, F[close], F_soft)
    Fvec = F[:, None] * rij
    N = x.shape[0]
    acc = np.zeros((N, 3))
    np.add.at(acc, ii,  Fvec * inv_mi[:, None])
    np.add.at(acc, jj, -Fvec * inv_mj[:, None])
    return acc

# ---------------------------------------------------------------------------
# Leapfrog integrator (general N). Fixed-step or heuristic adaptive sub-stepping.
# Returns times, pos, vel, info(force_evals, n_sub_history, ...)
# ---------------------------------------------------------------------------
def _prep(m, G):
    ii, jj = make_pairs(len(m))
    mi, mj = m[ii], m[jj]
    return dict(ii=ii, jj=jj, Gmimj=G*mi*mj, inv_mi=1.0/mi, inv_mj=1.0/mj)

def simulate_leapfrog(m, x0, v0, dt, T, n_samples, G, mode='newton',
                      adaptive=False, adapt_thresh=ADAPT_THRESH, max_substeps=MAX_SUBSTEPS):
    P = _prep(m, G); ii, jj = P['ii'], P['jj']
    eps2 = EPS*EPS
    x = x0.astype(np.float64).copy(); v = v0.astype(np.float64).copy()
    N = len(m)
    times = np.linspace(0.0, T, n_samples)
    n_steps = int(round(T/dt))
    pos = np.full((n_samples, N, 3), np.nan)   # NaN init: catch any unwritten slot
    vel = np.full((n_samples, N, 3), np.nan)
    fe = [0]
    def acc_of(xx):
        fe[0] += 1
        return compute_acc(xx, ii, jj, P['Gmimj'], P['inv_mi'], P['inv_mj'], eps2, mode)
    def min_sep(xx):
        d = xx[jj]-xx[ii]; return math.sqrt(float(np.min(np.einsum('ij,ij->i', d, d))) + 1e-30)
    a = acc_of(x)
    pos[0] = x; vel[0] = v
    si = 0
    n_sub_hist = []
    for step in range(n_steps):
        if adaptive and min_sep(x) < adapt_thresh:
            rmin = min_sep(x)
            n_sub = min(max_substeps, max(2, int(np.ceil(adapt_thresh/rmin))))
            h = dt/n_sub
            for _ in range(n_sub):
                vh = v + 0.5*h*a; x = x + h*vh; a = acc_of(x); v = vh + 0.5*h*a
            n_sub_hist.append(n_sub)
        else:
            vh = v + 0.5*dt*a; x = x + dt*vh; a = acc_of(x); v = vh + 0.5*dt*a
            n_sub_hist.append(1)
        t_cur = (step+1)*dt                       # index-derived: no accumulation drift
        while si < n_samples-1 and times[si+1] <= t_cur + 1e-9:
            si += 1; pos[si] = x; vel[si] = v
    while si < n_samples-1:                        # fill any trailing slot with final state
        si += 1; pos[si] = x; vel[si] = v
    # --- sampling guards ---
    assert not np.isnan(pos).any(), "unwritten sample slot remained (NaN)"
    assert not np.all(pos == 0.0, axis=2).all(axis=1).any(), "all-zero sample for all bodies"
    info = dict(force_evals=fe[0], n_sub_history=n_sub_hist,
                total_substeps=int(np.sum(n_sub_hist)), n_steps=n_steps,
                min_sep_run=float(min(min_sep(pos[k]) for k in range(0, n_samples, max(1, n_samples//200)))))
    return times, pos, vel, info

# ---------------------------------------------------------------------------
# REBOUND IAS15 reference
# ---------------------------------------------------------------------------
def ias15_reference(m, x0, v0, T, n_samples, G):
    import rebound
    sim = rebound.Simulation(); sim.integrator = "ias15"; sim.G = G
    for i in range(len(m)):
        sim.add(m=float(m[i]), x=float(x0[i,0]), y=float(x0[i,1]), z=float(x0[i,2]),
                vx=float(v0[i,0]), vy=float(v0[i,1]), vz=float(v0[i,2]))
    sim.move_to_com()
    times = np.linspace(0.0, T, n_samples)
    pos = np.zeros((n_samples, len(m), 3)); vel = np.zeros((n_samples, len(m), 3))
    for k, t in enumerate(times):
        sim.integrate(t)
        for i, p in enumerate(sim.particles):
            pos[k,i] = [p.x, p.y, p.z]; vel[k,i] = [p.vx, p.vy, p.vz]
    return times, pos, vel

# ---------------------------------------------------------------------------
# Diagnostics / metric panel
# ---------------------------------------------------------------------------
def energy_series(pos, vel, m, G, eps=EPS):
    KE = 0.5 * np.einsum("kij,i->k", vel**2, m)
    PE = np.zeros(pos.shape[0]); e2 = eps*eps; N = len(m)
    for i in range(N):
        for j in range(i+1, N):
            d = pos[:, i, :] - pos[:, j, :]; r2 = np.einsum("ki,ki->k", d, d)
            PE -= G * m[i]*m[j] / np.sqrt(r2 + e2)
    return KE + PE

def max_dE(pos, vel, m, G, eps=EPS):
    E = energy_series(pos, vel, m, G, eps)
    if not np.all(np.isfinite(E)): return float("inf")
    E0 = E[0]
    return float(np.max(np.abs((E - E0)/max(abs(E0), 1e-30))) * 100.0)

def rms_sep(a, b):
    d = a - b; pb = np.sqrt(np.sum(d**2, axis=-1))
    return np.sqrt(np.mean(pb**2, axis=1))

def bounded(pos, m):
    if not np.all(np.isfinite(pos)): return False
    com = (pos * m[None,:,None]).sum(1, keepdims=True) / m.sum()
    return bool(np.max(np.linalg.norm(pos - com, axis=2)) < EJECT_AU)

def metric_panel(times, pos, vel, pos_ref, m, G, force_evals, lyap=0.15, eps=EPS):
    """The panel of Table 3: energy, short-horizon accuracy, cost, bounded."""
    d = rms_sep(pos, pos_ref)
    T = times[-1]
    t_lyap = min(2.0/max(lyap, 1e-3), T)          # ~2 Lyapunov times
    short = d[times <= t_lyap]
    return dict(
        max_dE_pct = max_dE(pos, vel, m, G, eps),
        rms_short  = float(np.sqrt(np.mean(short**2))) if short.size else float("nan"),
        t_short    = float(t_lyap),
        rms_final  = float(d[-1]),                # coarse 100-yr sanity ONLY
        rms_timeavg= float(np.sqrt(np.mean(d**2))),
        force_evals= int(force_evals),
        bounded    = bounded(pos, m),
    )
