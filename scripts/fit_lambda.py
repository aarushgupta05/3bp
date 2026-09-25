"""fit_lambda.py - regenerate the four divergence rates of Table 1, unfloored.

The lambda values stored as literals in simon_core.py::_ICS_RAW were carried
over from earlier work and no script produced them. This is that script. It
reports the unfloored least-squares slope of log(separation) against time for
a shadow trajectory displaced by 1e-8, over the window 1e-7 < separation < 0.1,
with a 12-point robustness sweep, an exponential-versus-linear discrimination
and a renormalised Benettin cross-check.

Table 1 of the paper quotes the single-shadow fit column.

Run:  python scripts/fit_lambda.py [output_dir]
Out:  outputs/lambda_fit.txt, outputs/lambda_fit.json
"""
import os, sys, json, time, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)
import simon_core as sc

G = sc.G_TOY; T = 100.0; NS = 2001
CONFIGS = ["IC1", "IC3", "IC4", "IC6"]
D_LO, D_HI = 1e-7, 0.1
DELTAS = [1e-10, 1e-9, 1e-8, 1e-7]; BODIES = [0, 1, 2]
L = []
def log(s=""):
    print(s, flush=True); L.append(s)

def ias15(m,x0,v0,T,ns): return sc.ias15_reference(m,x0,v0,T,ns,G)
def fit_window(t,d): return (d>D_LO)&(d<D_HI)&np.isfinite(d)

def loglin_fit(t,d):
    msk=fit_window(t,d); n=int(msk.sum())
    if n<5: return dict(n=n,slope=float('nan'),se=float('nan'),r2=float('nan'),t0=float('nan'),t1=float('nan'))
    tt,y=t[msk],np.log(d[msk]); A=np.vstack([tt,np.ones_like(tt)]).T
    coef,*_=np.linalg.lstsq(A,y,rcond=None); pred=A@coef
    ss_res=float(np.sum((y-pred)**2)); ss_tot=float(np.sum((y-y.mean())**2))
    r2=1.0-ss_res/ss_tot if ss_tot>0 else float('nan')
    dof=max(n-2,1); s2=ss_res/dof; sxx=float(np.sum((tt-tt.mean())**2))
    se=math.sqrt(s2/sxx) if sxx>0 else float('nan')
    return dict(n=n,slope=float(coef[0]),se=float(se),r2=float(r2),t0=float(tt.min()),t1=float(tt.max()))

def lin_fit(t,d):
    msk=fit_window(t,d); n=int(msk.sum())
    if n<5: return dict(n=n,slope=float('nan'),r2=float('nan'))
    tt,y=t[msk],d[msk]; A=np.vstack([tt,np.ones_like(tt)]).T
    coef,*_=np.linalg.lstsq(A,y,rcond=None); pred=A@coef
    ss_res=float(np.sum((y-pred)**2)); ss_tot=float(np.sum((y-y.mean())**2))
    return dict(n=n,slope=float(coef[0]),r2=float(1.0-ss_res/ss_tot) if ss_tot>0 else float('nan'))

def shadow_separation(m,x0,v0,delta,body,ns=NS,horizon=T):
    t,p1,_=ias15(m,x0,v0,horizon,ns)
    xp=x0.copy(); xp[body,0]+=delta
    _,p2,_=ias15(m,xp,v0,horizon,ns)
    return t, sc.rms_sep(p1,p2)

def benettin(m,x0,v0,tau,delta=1e-8,horizon=T):
    import rebound
    def make(xx,vv):
        s=rebound.Simulation(); s.integrator='ias15'; s.G=G
        for i in range(len(m)):
            s.add(m=float(m[i]),x=float(xx[i,0]),y=float(xx[i,1]),z=float(xx[i,2]),
                  vx=float(vv[i,0]),vy=float(vv[i,1]),vz=float(vv[i,2]))
        return s
    def state(s):
        return (np.array([[p.x,p.y,p.z] for p in s.particles]),
                np.array([[p.vx,p.vy,p.vz] for p in s.particles]))
    s1=make(x0,v0); s1.move_to_com()
    x1,v1=state(s1)
    x2=x1.copy(); v2=v1.copy(); x2[0,0]+=delta
    s2=make(x2,v2)
    d0=delta
    n_steps=int(round(horizon/tau)); acc=0.0; used=0
    for k in range(n_steps):
        tgt=(k+1)*tau
        s1.integrate(tgt); s2.integrate(tgt)
        x1,v1=state(s1); x2,v2=state(s2)
        dx=x2-x1; dv=v2-v1
        d=math.sqrt(float(np.sum(dx**2)+np.sum(dv**2)))
        if d<=0 or not np.isfinite(d): break
        acc+=math.log(d/d0); used+=1
        f=d0/d
        s2=make(x1+dx*f, v1+dv*f); s2.t=s1.t
    return acc/(used*tau) if used else float('nan')

log("="*100)
log("REGENERATED DIVERGENCE RATES - pure Newtonian (point masses, no softening), G=1, T=100, NS=2001")
log("Unfloored least-squares slope of log(separation) vs t over 1e-7 < d < 0.1.")
log("="*100); log()
results={}; t_all=time.time()
for ic in CONFIGS:
    t0=time.time(); m,x0,v0=sc.get_ic(ic); stored=sc.ic_meta(ic)['exp_lambda']
    t,d=shadow_separation(m,x0,v0,1e-8,1)
    main=loglin_fit(t,d); lin=lin_fit(t,d)
    sweep=[]
    for dd in DELTAS:
        for b in BODIES:
            _,ds=shadow_separation(m,x0,v0,dd,b)
            sweep.append(dict(delta=dd,body=b,**loglin_fit(t,ds)))
    sl=np.array([s['slope'] for s in sweep],float); ok=sl[np.isfinite(sl)]
    ben={}
    for tau in (1.0,0.5): ben[f"tau={tau}"]=benettin(m,x0,v0,tau)
    results[ic]=dict(stored_literal=stored,main=main,linear=lin,sweep=sweep,
        sweep_median=float(np.median(ok)) if ok.size else float('nan'),
        sweep_min=float(ok.min()) if ok.size else float('nan'),
        sweep_max=float(ok.max()) if ok.size else float('nan'),
        benettin=ben,wall_s=time.time()-t0)
    r=results[ic]
    log(f"--- {ic} "+"-"*90)
    log(f"  stored literal in _ICS_RAW          : {stored:+.4f}")
    log(f"  single-shadow fit (d=1e-8, body 1)  : {main['slope']:+.6f}  +/- {main['se']:.2e}   R2(log-linear) = {main['r2']:.4f}   n = {main['n']}   window t = [{main['t0']:.2f}, {main['t1']:.2f}]")
    log(f"  linear fit of d(t), same window     : R2(linear)     = {lin['r2']:.4f}")
    better = "LOG-LINEAR (exponential)" if (np.isfinite(main['r2']) and np.isfinite(lin['r2']) and main['r2']>lin['r2']) else "LINEAR"
    log(f"  better model                        : {better}")
    log(f"  robustness sweep (12 estimates)     : median {r['sweep_median']:+.6f}   range [{r['sweep_min']:+.6f}, {r['sweep_max']:+.6f}]")
    for k,v in ben.items(): log(f"  Benettin renormalised, {k:<8}     : {v:+.6f}")
    log(f"  wall {r['wall_s']:.1f}s"); log()
log("="*100); log("VERDICT"); log("="*100)
for ic in CONFIGS:
    r=results[ic]; stored=r['stored_literal']; lo,hi=r['sweep_min'],r['sweep_max']
    inside=(lo<=stored<=hi) if np.isfinite(lo) and np.isfinite(hi) else False
    log(f"{ic}: stored {stored:+.4f} | sweep [{lo:+.6f}, {hi:+.6f}] | stored inside sweep range: {'YES' if inside else 'NO'}")
log()
r6=results['IC6']
log("IC6 SPECIFICALLY - is -0.0189 recovered?")
log(f"  single-shadow fit : {r6['main']['slope']:+.6f}")
log(f"  sweep range       : [{r6['sweep_min']:+.6f}, {r6['sweep_max']:+.6f}]")
log("  Benettin          : "+", ".join(f"{k} {v:+.6f}" for k,v in r6['benettin'].items()))
log(f"  R2 log-linear {r6['main']['r2']:.4f}  vs  R2 linear {r6['linear']['r2']:.4f}")
log(); log(f"total wall {time.time()-t_all:.1f}s")
open(os.path.join(OUT,"lambda_fit.txt"),"w").write("\n".join(L)+"\n")
json.dump(results,open(os.path.join(OUT,"lambda_fit.json"),"w"),indent=2,default=str)
print("\nwrote lambda_fit.txt + lambda_fit.json to",OUT)
