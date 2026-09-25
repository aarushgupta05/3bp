"""
integrators_symplectic.py - the time-symmetric adaptive leapfrog (TSALF) used
throughout the paper, alongside the fixed-step leapfrog in simon_core.py and
the REBOUND IAS15 reference.

TSALF follows Hut, Makino and McMillan (1995). The stepsize is solved
implicitly as h = 0.5 * (tau(start) + tau(end)) by a few fixed-point
iterations, which makes the step time-reversible, so the energy error stays
bounded even though the step varies. tau = eta * (minimum pairwise dynamical
time), where the dynamical time of a pair is the smaller of its orbital time
and its flyby time. eta is the one control parameter; the paper sweeps it over
0.2, 0.1, 0.05, 0.02, 0.01 and 0.005.

Energy is reported on the actual step states, not on the interpolated output
grid, so no smoothing hides a spike. Positions are interpolated onto the fixed
sample grid afterwards purely so the RMS comparison against IAS15 is made at
matching times.

Self-test. Run this file directly for a Kepler e = 0.9 energy-boundedness
check over 20 orbits and a forward-then-reversed round trip, which is the
property the implicit symmetric stepsize is supposed to buy.
"""
import math
import time as _time
import numpy as np
import simon_core as sc


def _acc_fn(m, G, mode):
    P = sc._prep(m, G); ii, jj = P['ii'], P['jj']; eps2 = sc.EPS * sc.EPS
    cnt = [0]
    def acc(x):
        cnt[0] += 1
        return sc.compute_acc(x, ii, jj, P['Gmimj'], P['inv_mi'], P['inv_mj'],
                              eps2, mode)
    return acc, cnt, P

def _verlet(x, v, a, h, acc):
    """One velocity-Verlet (KDK) step; 1 force eval. Returns x1,v1,a1."""
    vh = v + 0.5 * h * a
    x1 = x + h * vh
    a1 = acc(x1)
    v1 = vh + 0.5 * h * a1
    return x1, v1, a1

def t_dyn_min(x, v, m, ii, jj, G):
    """Min over pairs of min(orbital_time, flyby_time), the tightest timescale."""
    rij = x[jj] - x[ii]; vij = v[jj] - v[ii]
    r2 = np.einsum('ij,ij->i', rij, rij); r = np.sqrt(r2 + 1e-30)
    vmag = np.sqrt(np.einsum('ij,ij->i', vij, vij) + 1e-30)
    GM = G * (m[ii] + m[jj])
    t_orb = 2.0 * np.pi * np.sqrt(r2 * r / (GM + 1e-30))
    t_fly = r / vmag
    return float(np.min(np.minimum(t_orb, t_fly)))

# ---------------------------------------------------------------------------
def tsalf_simulate(m, x0, v0, eta, T, n_samples, G, mode='newton',
                   n_iter=2, max_steps=3_000_000, energy_eps=None, max_wall_s=None):
    acc, cnt, P = _acc_fn(m, G, mode); ii, jj = P['ii'], P['jj']
    x = x0.astype(np.float64).copy(); v = v0.astype(np.float64).copy(); a = acc(x)
    ts = [0.0]; xs = [x.copy()]; vs = [v.copy()]
    t = 0.0; nstep = 0
    # max_wall_s: abandon the run after this many seconds of wall time. The
    # partial result is still returned, with completed=False.
    _t_start = _time.time(); wall_stop = False
    while t < T and nstep < max_steps:
        if max_wall_s is not None and (nstep & 0xFFFF) == 0 and _time.time()-_t_start > max_wall_s:
            wall_stop = True; break
        tau0 = eta * t_dyn_min(x, v, m, ii, jj, G)
        h = tau0
        for _ in range(n_iter):                      # implicit symmetric stepsize
            x1, v1, a1 = _verlet(x, v, a, h, acc)
            tau1 = eta * t_dyn_min(x1, v1, m, ii, jj, G)
            h = 0.5 * (tau0 + tau1)
        if t + h > T:
            h = T - t                                # land exactly on T (last step only)
        x, v, a = _verlet(x, v, a, h, acc)
        t += h; nstep += 1
        ts.append(t); xs.append(x.copy()); vs.append(v.copy())
    ts = np.array(ts); xs = np.array(xs); vs = np.array(vs)
    # energy on the ACTUAL step states (no interpolation smoothing).
    # energy_eps=None keeps the softened default (eps=sc.EPS); pass
    # energy_eps=0.0 for a Newtonian run so the diagnostic matches the force model.
    maxdE = sc.max_dE(xs, vs, m, G, eps=(sc.EPS if energy_eps is None else energy_eps))
    # interpolate to the fixed sample grid for RMS-vs-IAS15
    times = np.linspace(0.0, T, n_samples)
    pos = np.empty((n_samples, len(m), 3)); vel = np.empty((n_samples, len(m), 3))
    for i in range(len(m)):
        for d in range(3):
            pos[:, i, d] = np.interp(times, ts, xs[:, i, d])
            vel[:, i, d] = np.interp(times, ts, vs[:, i, d])
    return times, pos, vel, dict(force_evals=cnt[0], n_steps=nstep, scheme='tsalf',
                                 eta=eta, max_dE_steps=maxdE, completed=(t >= T - 1e-9),
                                 wall_stopped=wall_stop, wall_s=_time.time()-_t_start,
                                 hit_step_cap=(nstep >= max_steps),
                                 step_states=(ts, xs, vs))

# ---------------------------------------------------------------------------
# Self-test: Kepler e=0.9, and the tsalf reversibility round trip.
# ---------------------------------------------------------------------------
def _kepler_ic(a=1.0, e=0.9, m0=1.0, m1=1e-3, G=1.0):
    M = m0 + m1
    rp = a * (1 - e); vp = math.sqrt(G * M * (1 + e) / rp)
    x = np.array([[0.0, 0, 0], [rp, 0, 0]]); v = np.array([[0.0, 0, 0], [0, vp, 0]])
    m = np.array([m0, m1])
    Mt = m.sum(); x -= (m[:, None] * x).sum(0) / Mt; v -= (m[:, None] * v).sum(0) / Mt
    P_orb = 2 * math.pi * math.sqrt(a**3 / (G * M))
    return m, x, v, P_orb

if __name__ == "__main__":
    G = 1.0
    m, x0, v0, Porb = _kepler_ic(a=1.0, e=0.9, G=G)
    T = 20 * Porb; NS = 4000
    print(f"=== KEPLER VERIFICATION (e=0.9, {20} orbits, Porb={Porb:.3f}) ===")
    for eta in [0.05, 0.02]:
        _, p, vv, info = tsalf_simulate(m, x0, v0, eta, T, NS, G, mode='newton')
        print(f"  tsalf    eta={eta:<5}: |dE/E0|max={info['max_dE_steps']:.3e}%  fe={info['force_evals']:7d}  steps={info['n_steps']}  bounded={sc.bounded(p,m)}")
    # reversibility: forward T then backward T (negate v), compare to start
    print("=== TSALF REVERSIBILITY (forward then time-reversed) ===")
    ts, xs, vs = tsalf_simulate(m, x0, v0, 0.05, 5*Porb, 2000, G, mode='newton')[3]['step_states']
    xf, vf = xs[-1].copy(), vs[-1].copy()
    _, _, _, infob = tsalf_simulate(m, xf, -vf, 0.05, 5*Porb, 2000, G, mode='newton')
    tb, xb, vb = infob['step_states']
    err = np.linalg.norm(xb[-1] - x0)
    print(f"  round-trip position error ||x_back - x_0|| = {err:.3e}  (small => reversible)")
    print("KEPLER_SELFTEST_DONE")
