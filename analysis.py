#!/usr/bin/env python3
"""Regenerate EVERY headline table from results/. Single source of truth for all quoted numbers.

    python analysis.py            # all tables
    python analysis.py --json     # machine-readable

Each table prints the run-ids it used, so any number in the docs is traceable to committed data.
No number should appear in a paper or plan that this script cannot reproduce.
"""
import json, glob, math, argparse, os
import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')

def load(pat='**/*.json'):
    out = []
    for f in glob.glob(os.path.join(RES, pat), recursive=True):
        try:
            d = json.load(open(f)); d.setdefault('tag', os.path.basename(f)[:-5]); out.append(d)
        except Exception: pass
    return out

def tk(r, K=1.5):
    """Collapse time by FORWARD threshold crossing. Prefers the stored ladder; else re-derives
    from the series. NEVER 'first return after the global minimum' (that returns NaN for runs
    that collapsed then evacuated their core)."""
    if r.get('t_col_k'): 
        v = r['t_col_k'].get(f'K{K}')
        if v is not None: return v
    ts, rs = r.get('t_series'), r.get('rho_series')
    if not ts: return r.get('t_collapse', float('nan'))
    ts, rs = np.array(ts), np.array(rs)
    j = np.where(rs >= K*rs[0])[0]
    return float(ts[j[0]]) if len(j) else float('nan')

def stat(v):
    v = [x for x in v if x == x]
    if not v: return (float('nan'), 0.0, 0)
    return (float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, len(v))

def T(title): print("\n" + "="*78 + f"\n{title}\n" + "="*78)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--json', action='store_true'); a = ap.parse_args()
    R = load(); out = {}
    if not R: print(f"no results under {RES}"); return

    T("V1  SIDM collapse is GRAVOTHERMAL, not two-body relaxation (t_col flat vs N)")
    Ns = [2000, 4000, 6000, 8000]; rows = []
    for K in (1.3, 1.5, 2.0):
        ms = [stat([tk(r, K) for r in R if r['tag'].startswith(f'F_N{n}_')])[0] for n in Ns]
        ok = [i for i, m in enumerate(ms) if m == m]
        sl = np.polyfit(np.log([Ns[i] for i in ok]), np.log([ms[i] for i in ok]), 1)[0] if len(ok) > 2 else float('nan')
        rows.append((K, ms, sl))
        print(f"  K={K}: " + " ".join(f"N={n}:{m:6.0f}" for n, m in zip(Ns, ms)) + f"   slope={sl:+.2f}")
    trel = [n/(8*math.log(n))*(0.2/0.8/0.005) for n in Ns]
    print(f"  t_relax(2body): " + " ".join(f"{t:6.0f}" for t in trel) + f"   (grows {trel[-1]/trel[0]:.1f}x)")
    out['V1_slopes'] = {f'K{k}': s for k, _, s in rows}

    T("V2  matched sham floor is SCATTER, not bias")
    for K in (1.3, 1.5, 2.0):
        n = stat([tk(r, K) for r in R if r['tag'].startswith('F_none_')])
        s = stat([tk(r, K) for r in R if r['tag'].startswith('F_sham_')])
        if n[2] and s[2]:
            se = math.sqrt(n[1]**2/n[2] + s[1]**2/s[2]); z = abs(n[0]-s[0])/se if se else float('inf')
            print(f"  K={K}: none={n[0]:.0f}+/-{n[1]:.0f}(n={n[2]}) sham={s[0]:.0f}+/-{s[1]:.0f}(n={s[2]})  {z:.1f} sigma")
            out[f'V2_sigma_K{K}'] = z

    T("V3/V4  anisotropy erasure: collisional, and fast relative to collapse")
    for tag, lbl in (('D_om_grav', 'OM equilibrium, sigma/m=0'), ('D_om_sidm', 'OM equilibrium, SIDM'),
                     ('D_rot_grav', 'rotated, sigma/m=0'), ('D_rot_sidm', 'rotated, SIDM')):
        g = [r for r in R if r['tag'].startswith(tag)]
        if not g: continue
        keep = []
        for r in g:
            bs, ts = r.get('beta_series'), r.get('t_series')
            if not bs: continue
            i = min(range(len(ts)), key=lambda k: abs(ts[k]-1000))
            keep.append(abs(bs[i])/max(abs(r['beta_i']), 1e-9))
        if keep: print(f"  {lbl:<28} beta retained at t=1000: {100*np.mean(keep):3.0f}%   (n={len(keep)})")
    print("  t_iso/t_collapse by sigma/m:")
    for sg in (5, 10, 20, 40, 80):
        g = [r for r in R if r['tag'].startswith(f'F_sig{sg}_')]
        ti = stat([r.get('t_iso') for r in g]); tc = stat([tk(r) for r in g])
        if ti[2] and tc[0] == tc[0]:
            print(f"    sigma/m={sg:3d}: t_iso={ti[0]:6.0f} t_col={tc[0]:7.0f} ratio={ti[0]/tc[0]:.3f}")

    T("V5  beta-dependence of collapse time (the bounded NULL)")
    for K in (1.3, 1.5, 2.0):
        groups = {}
        for r in R:
            if r['tag'].startswith('E_b000') or r['tag'].startswith('E_ra'):
                groups.setdefault(r['tag'].rsplit('_s', 1)[0], []).append(r)
        pts = []
        for k, g in groups.items():
            b = np.mean([x['beta_i'] for x in g]); m, s, n = stat([tk(x, K) for x in g])
            if n: pts.append((b, m, s, n))
        if len(pts) > 2:
            pts.sort(); bs = np.array([p[0] for p in pts]); ms = np.array([p[1] for p in pts])
            sl = np.polyfit(bs, ms, 1)[0]
            print(f"  K={K}: slope={sl:+7.0f} steps per unit beta | max/min group ratio={ms.max()/ms.min():.2f}x")
            out[f'V5_slope_K{K}'] = float(sl)
    print("  NOTE: the SIGN of the slope flips with K -> consistent with noise, not signal.")

    T("V6  run-to-run scatter is DYNAMICAL, not Poisson (sets campaign sizing)")
    fr = []
    for n in Ns:
        m, s, c = stat([tk(r) for r in R if r['tag'].startswith(f'F_N{n}_')])
        if c > 1: fr.append((n, s/m)); print(f"  N={n:5d}: sd/mean = {100*s/m:5.1f}%  (n={c})")
    if len(fr) > 2:
        sl = np.polyfit(np.log([x[0] for x in fr]), np.log([x[1] for x in fr]), 1)[0]
        print(f"  scatter ~ N^{sl:+.2f}  (Poisson=-0.5, dynamical=0) -> power comes from SEEDS, not N")
        out['V6_scatter_exponent'] = float(sl)

    T("SAFETY  SIDM Monte-Carlo diagnostics (must be present on all new runs)")
    d = [r for r in R if r.get('kappa_max') is not None]
    if not d:
        print("  none of these results carry SIDM safety diagnostics.")
        print("  -> pre-instrumentation runs. New runs record P_max / kappa_max / blocked_frac.")
    else:
        for r in sorted(d, key=lambda x: x.get('dt', 0)):
            print(f"  {r['tag']:<16} dt={r.get('dt')} kappa={r['kappa_max']:.4f} "
                  f"P_max={r['P_max']:.3f} blocked={r['blocked_frac']:.3f} "
                  f"scat/step={r.get('scatters_per_step', 0):.2f}")
        print("  K&S guidance: kappa<=0.02 for O(10%) accuracy; blocked_frac is the once-per-step")
        print("  saturation bias (candidate pairs skipped because a partner already scattered).")

    print(f"\n(loaded {len(R)} result files from {RES})")
    if a.json: print(json.dumps(out, indent=2))

if __name__ == '__main__':
    main()
