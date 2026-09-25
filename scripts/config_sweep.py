"""
config_sweep.py - the configuration sweep behind Table 2 and Fig. 2.

Ten configurations (IC1/IC3/IC4/IC6, binary-single at r_bin 0.05/0.10/0.20/0.40,
and Sun-Earth-Moon from JPL Horizons) against six method variants (heuristic
adaptive leapfrog, fixed-step leapfrog at dt = 0.04/0.02/0.01/0.005, TSALF at
eta = 0.05, and IAS15). Metrics recorded per run: force evals per year, wall
per year, max |dE/E0|, bounded, and the pair-separation range where a pair is
expected to stay bound.

Force model: pure Newtonian by default, which is what the paper uses, with the
energy diagnostic at eps = 0 so the metric matches the force model. Set
SIMON_FORCE_MODEL=soft for the softened model; that writes to a separate
filename so a Newtonian output is never overwritten.

Sun-Earth-Moon is skipped by default. It is already Newtonian in either force
model, none of Fig. 2 or Table 2 needs its rows, and running it would refetch
live state vectors from JPL Horizons, which would put the paper's Fig. 4
headline at the mercy of a network call. The frozen Sun-Earth-Moon run that
Fig. 4 uses is committed at outputs/sun_earth_moon_run.json. Set SIMON_SKIP_SEM=0
to refetch and rerun it, and expect the numbers to move if Horizons has
changed its ephemeris.

Run:  python scripts/config_sweep.py
Out:  outputs/config_sweep_newton.json
"""
import sys, os, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import numpy as np
import simon_core as sc
import integrators_symplectic as si

OUT = os.path.join(ROOT, "outputs"); os.makedirs(OUT, exist_ok=True)
T, NS = 100.0, 5000
TSALF_CAP = 2_500_000

# --- force model ----------------------------------------------------------
# 'newton' (default) is what the paper uses: pure Newtonian for IC1/IC3/IC4/IC6
# and every BS_* config, with the ENERGY DIAGNOSTIC at eps=0 so the metric
# matches the force model. 'soft' keeps the softened model available. The
# Sun-Earth-Moon configuration, which the sweep labels SEM, is Newtonian
# either way.
FORCE_MODEL = os.environ.get("SIMON_FORCE_MODEL", "newton")
assert FORCE_MODEL in ("soft", "newton"), FORCE_MODEL
E_EPS = 0.0 if FORCE_MODEL == "newton" else sc.EPS
SUFFIX = "_newton" if FORCE_MODEL == "newton" else "_soft"
# Sun-Earth-Moon is skipped by default: it is already Newtonian, neither
# Fig. 2 nor Table 2 needs its rows, and skipping avoids a live JPL Horizons
# refetch, which would be a reproducibility risk for no benefit. The frozen
# run is at outputs/sun_earth_moon_run.json.
SKIP_SEM = os.environ.get("SIMON_SKIP_SEM", "1") == "1"

# leapfrog_0.02 covers the band where fixed-step leapfrog ejects, which the
# paper's Finding 1 depends on.
METHODS = ["heuristic", "leapfrog_0.04", "leapfrog_0.02", "leapfrog_0.01",
           "leapfrog_0.005", "tsalf", "ias15"]

def binary_single_ic(r_bin, G):
    m = np.array([1.0, 0.5, 0.1]); mb = 1.5
    x = np.array([[-r_bin*0.5/mb,0,0],[r_bin*1.0/mb,0,0],[3.0,0,0]], float)
    vrel = math.sqrt(G*mb/r_bin); v2 = math.sqrt(G*mb/3.0)
    v = np.array([[0,-vrel*0.5/mb,0],[0,vrel*1.0/mb,0],[0,v2,0]], float)
    M = m.sum(); x -= (m[:,None]*x).sum(0)/M; v -= (m[:,None]*v).sum(0)/M
    return m, x, v

def sem_ic():
    from astroquery.jplhorizons import Horizons
    def vec(bid):
        o = Horizons(id=bid, location="@0", epochs={"start":"2026-01-01","stop":"2026-01-02","step":"1d"})
        t = o.vectors()
        p = np.array([float(t["x"][0]), float(t["y"][0]), float(t["z"][0])])
        vv = np.array([float(t["vx"][0]), float(t["vy"][0]), float(t["vz"][0])]) * 365.25
        return p, vv
    xs,vs = vec("10"); xe,ve = vec("399"); xm,vm = vec("301")
    m = np.array([1.0, 3.003489614915e-6, 3.694303349e-8])
    x0 = np.vstack([xs,xe,xm]); v0 = np.vstack([vs,ve,vm])
    M = m.sum(); x0 -= (m[:,None]*x0).sum(0)/M; v0 -= (m[:,None]*v0).sum(0)/M
    return m, x0, v0

def ias15_run(m, x0, v0, G):
    import rebound
    sim = rebound.Simulation(); sim.integrator = "ias15"; sim.G = G
    for i in range(len(m)):
        sim.add(m=float(m[i]), x=float(x0[i,0]), y=float(x0[i,1]), z=float(x0[i,2]),
                vx=float(v0[i,0]), vy=float(v0[i,1]), vz=float(v0[i,2]))
    sim.move_to_com()
    times = np.linspace(0, T, NS); pos = np.zeros((NS,len(m),3)); vel = np.zeros((NS,len(m),3))
    for k,t in enumerate(times):
        sim.integrate(t)
        for i,p in enumerate(sim.particles):
            pos[k,i]=[p.x,p.y,p.z]; vel[k,i]=[p.vx,p.vy,p.vz]
    steps = int(sim.steps_done)
    return pos, vel, steps, steps*8   # fe estimate: 8 Gauss-Radau nodes/step (labeled)

def pair_range(pos, idx):
    d = np.linalg.norm(pos[:,idx[0],:]-pos[:,idx[1],:], axis=1); return float(d.min()), float(d.max())

# ---- config registry ----
CONFIGS = []
for ic in ["IC1","IC3","IC4","IC6"]:
    m,x0,v0 = sc.get_ic(ic)
    CONFIGS.append(dict(name=ic, cls="existing", G=sc.G_TOY, mode=FORCE_MODEL, op_dt=0.04, pair=None, m=m, x0=x0, v0=v0))
for rb in [0.05,0.10,0.20,0.40]:
    m,x0,v0 = binary_single_ic(rb, sc.G_REAL)
    CONFIGS.append(dict(name=f"BS_{rb:.2f}", cls="binary", G=sc.G_REAL, mode=FORCE_MODEL, op_dt=0.04, pair=(0,1), r_bin=rb, m=m, x0=x0, v0=v0))
if SKIP_SEM:
    print("[sweep] Sun-Earth-Moon SKIPPED (SIMON_SKIP_SEM=1): already Newtonian, no "
          "Horizons refetch. Fig. 4 uses outputs/sun_earth_moon_run.json.", flush=True)
else:
    print("[sweep] fetching Sun-Earth-Moon initial conditions from JPL Horizons ...", flush=True)
    m,x0,v0 = sem_ic()
    CONFIGS.append(dict(name="SEM", cls="sem", G=sc.G_REAL, mode="newton", op_dt=0.016, pair=(1,2), m=m, x0=x0, v0=v0))

results = []; sem_series = {}
for cfg in CONFIGS:
    m,x0,v0,G,mode = cfg["m"],cfg["x0"],cfg["v0"],cfg["G"],cfg["mode"]
    print(f"[sweep] {cfg['name']} ({cfg['cls']}) ...", flush=True)
    for meth in METHODS:
        t0 = time.time()
        try:
            if meth == "heuristic":
                _,p,v,info = sc.simulate_leapfrog(m,x0,v0,cfg["op_dt"],T,NS,G,mode=mode,adaptive=True)
                fe = info["force_evals"]; done = T; maxdE = sc.max_dE(p,v,m,G,eps=E_EPS)
            elif meth.startswith("leapfrog_"):
                dt = float(meth.split("_")[1])
                _,p,v,info = sc.simulate_leapfrog(m,x0,v0,dt,T,NS,G,mode=mode,adaptive=False)
                fe = info["force_evals"]; done = T; maxdE = sc.max_dE(p,v,m,G,eps=E_EPS)
            elif meth == "tsalf":
                _,p,v,info = si.tsalf_simulate(m,x0,v0,0.05,T,NS,G,mode=mode,
                                               max_steps=TSALF_CAP,energy_eps=E_EPS)
                fe = info["force_evals"]; done = T if info["completed"] else float(info["step_states"][0][-1])
                maxdE = info["max_dE_steps"]
            elif meth == "ias15":
                p,v,steps,fe = ias15_run(m,x0,v0,G); done = T; maxdE = sc.max_dE(p,v,m,G,eps=E_EPS)
            wall = time.time()-t0
            bnd = sc.bounded(p,m)
            rec = dict(config=cfg["name"], cls=cfg["cls"], method=meth,
                       fe_per_yr=fe/max(done,1e-9), wall_per_yr=wall/max(done,1e-9),
                       max_dE_pct=(maxdE if np.isfinite(maxdE) else None), bounded=bnd,
                       completed_yr=round(done,3), total_fe=int(fe), wall_s=round(wall,3))
            if cfg["pair"] is not None:
                pr = pair_range(p, cfg["pair"]); rec["pair_range"] = [pr[0], pr[1]]
            if cfg["cls"] == "sem":
                em = np.linalg.norm(p[:,1,:]-p[:,2,:], axis=1)
                sem_series[meth] = dict(t=np.linspace(0,T,NS)[::5].tolist(), em=em[::5].tolist())
            results.append(rec)
            dEs = "FAIL" if rec["max_dE_pct"] is None else f"{rec['max_dE_pct']:.3g}%"
            print(f"   {meth:<14} fe/yr={rec['fe_per_yr']:.1f} dE={dEs} bnd={bnd} done@{done:.1f}yr wall={wall:.1f}s", flush=True)
        except Exception as e:
            results.append(dict(config=cfg["name"], cls=cfg["cls"], method=meth, error=str(e)))
            print(f"   {meth:<14} ERROR {e}", flush=True)

payload = dict(T=T, NS=NS, force_model=FORCE_MODEL, energy_eps=E_EPS,
               skipped_sem=SKIP_SEM, methods=METHODS,
               configs=[{k:c[k] for k in ("name","cls","G","mode","op_dt") } for c in CONFIGS],
               results=results, sem_series=sem_series)
# ground rule: the two force models never overwrite each other
dest = os.path.join(OUT, f"config_sweep{SUFFIX}.json")
json.dump(payload, open(dest,"w"), indent=2, default=str)
print(f"[sweep] wrote {dest}", flush=True)
print("CONFIG_SWEEP_DONE", flush=True)
