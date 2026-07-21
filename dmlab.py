#!/usr/bin/env python3
"""
dmlab -- ONE dark-matter experiment lab: shared physics core + experiment registry.

Replaces the scattered per-idea scripts (observer_gamma, peri_handle, collapse_test, dm_residual,
dm_scaleflow, causal_locus, shape_observer, loschmidt*, mond_response, superfluid_response,
kerr_losscone, stream_shape), which had `gamma` duplicated across 8 files and `radialize` across 8.

  python dmlab.py list                      # registered experiments
  python dmlab.py run <name> [--N ..] [--steps ..] [--seeds ..]
  python dmlab.py verify                    # regression: reproduce recorded published numbers

Core sections: SAMPLING | INTERVENTIONS | MEASUREMENT | FORCE LAWS + EVOLVE | EXPERIMENTS | CLI.
Units G = M_halo = 1 unless noted. Builds on validated modules: nbody_3d (acceleration),
nbody_dm_ic (sample_nfw3d), nfw_anisotropic_ic (equilibrium OM ICs), sidm, feedback.
"""
from __future__ import annotations
import argparse, math, types, sys
import numpy as np

from nbody_3d import SimConfig, acceleration
from nbody_dm_ic import sample_nfw3d
from nfw_anisotropic_ic import sample_nfw_anisotropic

CEN = np.full(3, 1.0)          # halo centre (box_size=2.0)
KAPPA_TARGET = 0.0143          # max scatters/particle per SIDM sub-step.
# NOTE: sidm_scatter's recorded kappa UNDERCOUNTS by a factor 0.716 (it accumulates P only
# over pairs it evaluates, skipping blocked ones). The K&S criterion is <=0.02 on the TRUE
# value, so the target applied here must be 0.02*0.716 = 0.0143. With the old 0.02 target the
# true kappa landed at 0.02/0.716 = 0.0279 -- confirmed at 0.0270-0.0273 across every conv5
# cell, i.e. `auto` was systematically under-sub-cycling by 40%.
RS_DEF, C_DEF, EPS_DEF, DT_DEF = 0.20, 10.0, 0.05, 0.005


# ────────────────────────────── SAMPLING ──────────────────────────────

def cfg_for(N, eps=EPS_DEF, rs=RS_DEF, dt=DT_DEF):
    return SimConfig(model='direct_isolated', integrator='leapfrog_kdk', init='nfw3d',
                     seed=1, n=N, eps=eps, box_size=2.0, plummer_a=rs, steps=1, dt=dt)

def sample(N, seed=1, rs=RS_DEF, c=C_DEF, kind='isotropic', r_a=None, eps=EPS_DEF):
    """kind: 'nfw3d' (fast isotropic sampler) | 'isotropic'/'radial' (equilibrium OM Jeans ICs)."""
    if kind == 'nfw3d':
        icfg = types.SimpleNamespace(n=N, plummer_a=rs, nfw_c=c, box_size=2.0)
        return sample_nfw3d(np.random.default_rng(seed), icfg)
    return sample_nfw_anisotropic(N, rs, c, kind, np.random.default_rng(seed),
                                  r_a=r_a, box_size=2.0, eps=eps)


# ─────────────────────────── INTERVENTIONS ────────────────────────────
# All speed-preserving (positions, per-particle speed, KE exactly fixed) and momentum re-zeroed.
# `mask` restricts the intervention to a particle subset (e.g. a radial shell).

def _frame(pos, vel):
    d = pos - CEN; r = np.linalg.norm(d, axis=1)
    rhat = d / np.maximum(r, 1e-12)[:, None]
    v_r = np.sum(vel * rhat, axis=1)
    vt = vel - v_r[:, None] * rhat
    v_t = np.linalg.norm(vt, axis=1)
    return rhat, v_r, vt, v_t, np.linalg.norm(vel, axis=1)

def _apply(vel, vnew, mask):
    out = vel.copy()
    out[mask] = vnew[mask]
    return out - np.mean(out, axis=0)

def radialize(pos, vel, theta, mask=None):
    """Rotate velocities toward the NEAREST radial direction by theta (raises f_peri, beta up)."""
    rhat, v_r, vt, v_t, sp = _frame(pos, vel)
    that = vt / np.where(v_t > 1e-12, v_t, 1.0)[:, None]
    phi = np.arctan2(v_t, v_r)
    phi2 = np.where(phi <= math.pi/2, np.maximum(phi - theta, 0.0), np.minimum(phi + theta, math.pi))
    vnew = sp[:, None]*(np.cos(phi2)[:, None]*rhat + np.sin(phi2)[:, None]*that)
    m = (np.ones(len(vel), bool) if mask is None else mask) & (v_t > 1e-9)
    return _apply(vel, vnew, m)

def tangentialize(pos, vel, theta, mask=None, seed=3):
    """Rotate velocities TOWARD the tangential plane by theta (lowers f_peri, beta down)."""
    rhat, v_r, vt, v_t, sp = _frame(pos, vel)
    that = np.zeros_like(vel); big = v_t > 1e-9
    that[big] = vt[big] / v_t[big][:, None]
    if np.any(~big):                                   # purely radial: pick a random t-hat perp r-hat
        rng = np.random.default_rng(seed ^ 0x7777)
        rnd = rng.normal(size=(int((~big).sum()), 3)); rs_ = rhat[~big]
        rnd -= np.sum(rnd*rs_, axis=1)[:, None]*rs_
        rnd /= np.linalg.norm(rnd, axis=1)[:, None]; that[~big] = rnd
    phi = np.arctan2(v_t, v_r)
    phi2 = np.where(phi <= math.pi/2, np.minimum(phi + theta, math.pi/2), np.maximum(phi - theta, math.pi/2))
    vnew = sp[:, None]*(np.cos(phi2)[:, None]*rhat + np.sin(phi2)[:, None]*that)
    return _apply(vel, vnew, np.ones(len(vel), bool) if mask is None else mask)

def sham(pos, vel, theta, mask=None, seed=11):
    """Random-axis rotation: speed-preserving and anisotropy-NEUTRAL. The matched null control."""
    rng = np.random.default_rng(seed)
    k = rng.normal(size=vel.shape); k /= np.linalg.norm(k, axis=1, keepdims=True)
    kxv = np.cross(k, vel); kd = np.sum(k*vel, axis=1)
    vnew = vel*math.cos(theta) + kxv*math.sin(theta) + k*kd[:, None]*(1-math.cos(theta))
    return _apply(vel, vnew, np.ones(len(vel), bool) if mask is None else mask)

def l_null(pos, vel, phi_rot, mask=None):
    """Reorient the ORBITAL PLANE: rotate v_t about r-hat by phi_rot. Preserves speed, |L|, E and
    therefore PERICENTER exactly -- changes only L-hat. The instrument that is null for spherical
    systems but causal for orientation-dependent (Kerr) boundaries."""
    rhat, v_r, vt, v_t, sp = _frame(pos, vel)
    that = np.zeros_like(vel); big = v_t > 1e-9
    that[big] = vt[big] / v_t[big][:, None]
    perp = np.cross(rhat, that)                        # completes (rhat, that, perp)
    vnew = v_r[:, None]*rhat + v_t[:, None]*(math.cos(phi_rot)*that + math.sin(phi_rot)*perp)
    m = (np.ones(len(vel), bool) if mask is None else mask) & big
    return _apply(vel, vnew, m)


# ──────────────────────────── MEASUREMENT ─────────────────────────────

def gamma(pos, lo=0.05, hi=0.35, nb=7, minc=12):
    """Inner logarithmic density slope -dln(rho)/dln(r) over [lo,hi]."""
    r = np.linalg.norm(pos - CEN, axis=1)
    e = np.logspace(np.log10(lo), np.log10(hi), nb+1); rc, rho = [], []
    for i in range(nb):
        m = ((r >= e[i]) & (r < e[i+1])).sum()
        v = 4/3*np.pi*(e[i+1]**3 - e[i]**3)
        if m > minc: rc.append((e[i]*e[i+1])**.5); rho.append(m/v)
    return -np.polyfit(np.log(rc), np.log(rho), 1)[0]

def pericenters(pos, vel):
    """Per-particle pericenter in the halo's OWN spherically-averaged potential (self-consistent)."""
    N = len(pos)
    d = pos - CEN; r = np.linalg.norm(d, axis=1); v2 = np.sum(vel*vel, axis=1)
    L = np.linalg.norm(np.cross(d, vel), axis=1)
    rg = np.logspace(np.log10(max(r.min(), 1e-3)), np.log10(r.max()), 500)
    Mr = np.searchsorted(np.sort(r), rg)/N
    integ = Mr/rg**2
    Phi = -(Mr[-1]/rg[-1]) - np.concatenate(
        [np.cumsum((0.5*(integ[1:]+integ[:-1])*np.diff(rg))[::-1])[::-1], [0.0]])
    E = 0.5*v2 + np.interp(r, rg, Phi)
    g = (Phi[None, :] + (L[:, None]**2)/(2.0*rg[None, :]**2)) - E[:, None]
    below = g < 0; rp = np.full(N, rg[0]); has = below.any(axis=1)
    rp[has] = rg[np.argmax(below, axis=1)[has]]
    return rp

def f_peri(pos, vel, r_c=0.10):
    return float(np.mean(pericenters(pos, vel) < r_c))

def circularity(pos, vel, rs=RS_DEF, c=C_DEF):
    """eta = L / L_circ(E) in the truncated-NFW potential. Radius-INDEPENDENT orbit-type measure
    (a raw-L split is confounded: small-r particles mechanically have small L)."""
    d = pos - CEN; r = np.linalg.norm(d, axis=1); v2 = np.sum(vel*vel, axis=1)
    L = np.linalg.norm(np.cross(d, vel), axis=1)
    def g(x): return np.log1p(x) - x/(1.0+x)
    gc = g(c)
    Mlt = lambda rr: g(np.clip(rr/rs, 0, c))/gc
    def Phi(rr):
        x = rr/rs; return -(1.0/(gc*rs))*(np.log1p(x)/x - 1.0/(1.0+c))
    rg = np.logspace(np.log10(1e-3*rs), np.log10(c*rs*0.999), 4000)
    vc2 = Mlt(rg)/rg; Ecirc = 0.5*vc2 + Phi(rg); Lcirc = rg*np.sqrt(vc2)
    E = np.clip(0.5*v2 + Phi(r), Ecirc[0], Ecirc[-1])
    Lc = np.interp(np.interp(E, Ecirc, rg), rg, Lcirc)
    return np.clip(L/np.maximum(Lc, 1e-30), 0, 1.5)

def beta_sigma_profiles(pos, vel, edges):
    """Binned rho(r), sigma_r^2(r), anisotropy beta(r), and true M(<r) from particle counts."""
    N = len(pos); d = pos - CEN; r = np.linalg.norm(d, axis=1)
    rhat = d/np.maximum(r, 1e-12)[:, None]
    vr = np.sum(vel*rhat, axis=1); v2 = np.sum(vel*vel, axis=1)
    rc = np.sqrt(edges[:-1]*edges[1:]); shellV = 4/3*np.pi*(edges[1:]**3 - edges[:-1]**3)
    rho = np.empty(len(rc)); sr2 = np.empty(len(rc)); bet = np.empty(len(rc))
    for i in range(len(rc)):
        m = (r >= edges[i]) & (r < edges[i+1]); n = m.sum()
        rho[i] = n/N/shellV[i]
        vri = vr[m]; sr2[i] = np.mean(vri**2) - np.mean(vri)**2
        bet[i] = 1.0 - (np.mean(v2[m]) - np.mean(vri**2))/(2.0*sr2[i])
    return rc, rho, sr2, bet, np.searchsorted(np.sort(r), rc)/N

BETA_EDGES = np.logspace(np.log10(0.05), np.log10(0.80), 16)
BETA_SEL = slice(2, 10)

def mean_beta(pos, vel):
    return float(np.mean(beta_sigma_profiles(pos, vel, BETA_EDGES)[3][BETA_SEL]))

def dose_for_beta(pos, vel, beta_target, tol=0.02, itmax=40):
    """Solve the rotation angle that yields a TARGET mean anisotropy beta.

    Necessary because the dose->beta map is strongly ASYMMETRIC: beta is bounded above by 1 but
    unbounded below (sigma_r -> 0), so theta=0.8 gives beta=+0.92 radial but beta=-17 tangential.
    Specifying theta directly would produce mismatched arms and a pathological tangential state.
    Returns (arm, theta). Bisection on the monotonic branch.
    """
    if abs(beta_target) < 1e-6: return ('none', 0.0)
    arm = 'radial' if beta_target > 0 else 'tang'
    rot = (lambda th: radialize(pos, vel, th)) if arm == 'radial' else (lambda th: tangentialize(pos, vel, th))
    lo, hi = 0.0, 1.4
    for _ in range(itmax):
        mid = 0.5*(lo+hi); b = mean_beta(pos, rot(mid))
        if abs(b - beta_target) < tol: return (arm, mid)
        # beta increases with theta for radial, decreases for tangential
        if (b < beta_target) == (arm == 'radial'): lo = mid
        else: hi = mid
    return (arm, 0.5*(lo+hi))

def jeans_mass(rc, rho, sr2, beta):
    """Spherical Jeans M(<r). Returns (M_isotropic_assumed, M_beta_aware).
    A beta-blind observer drops the 2*beta term -> that difference IS the observer-relative residual."""
    dlnP = np.gradient(np.log(rho*sr2), np.log(rc))
    return -(rc*sr2)*dlnP, -(rc*sr2)*(dlnP + 2.0*beta)

def axis_ratios(pts, rmax=None):
    """(c/a, b/a) from the inertia tensor. NOTE: for near-planar orbit footprints c/a is dominated by
    ORBITAL PLANARITY, not halo shape -- use b/a (the in-plane ratio) as the shape estimator there."""
    q = pts - np.median(pts, axis=0)
    if rmax is not None: q = q[np.linalg.norm(q, axis=1) < rmax]
    w = np.sort(np.linalg.eigvalsh(q.T @ q / len(q)))[::-1]
    return math.sqrt(w[2]/w[0]), math.sqrt(w[1]/w[0])


# ──────────────────────── FORCE LAWS + EVOLVE ─────────────────────────

def make_force(law, cfg, mass, a0=1.0, r_sf=0.5):
    """law: 'newton' | 'mond' (algebraic Milgrom nu) | 'superfluid' (MOND gated to r<r_sf)."""
    def newton(p): return acceleration(p, mass, cfg, use_numba=True)
    if law == 'newton': return newton
    def nu(aN):
        gN = np.linalg.norm(aN, axis=1, keepdims=True)
        return 0.5 + np.sqrt(0.25 + a0/np.maximum(gN, 1e-12))
    if law == 'mond':
        return lambda p: nu(newton(p))*newton(p)
    if law == 'superfluid':
        def sf(p):
            aN = newton(p); inside = (np.linalg.norm(p - CEN, axis=1) < r_sf)[:, None]
            return np.where(inside, nu(aN)*aN, aN)
        return sf
    raise ValueError(law)

def evolve(pos, vel, cfg, mass, steps, law='newton', dt=DT_DEF, a0=1.0, r_sf=0.5,
           sigma_over_m=0.0, h=0.10, rng=None, feedback_amp=0.0, r_gas=0.10, T_osc=120,
           sidm_diag=None, sidm_subcycles='auto'):
    """Unified leapfrog KDK. Optional SIDM scattering and/or oscillating-central-gas feedback."""
    accel = make_force(law, cfg, mass, a0, r_sf)
    from sidm import sidm_scatter
    rng = rng or np.random.default_rng(7)
    p, v = pos.copy(), vel.copy(); a = accel(p); n_scat = 0
    if sidm_diag is None: sidm_diag = {}
    sidm_diag.setdefault('n_scatter_total', 0)
    for s in range(steps):
        if feedback_amp:
            Mg = feedback_amp*(0.5 - 0.5*math.cos(2*math.pi*s/T_osc))
            dx = p - CEN; rr2 = np.sum(dx*dx, axis=1)[:, None] + r_gas**2
            a = a - Mg*dx/rr2**1.5
        v = v + 0.5*dt*a; p = p + dt*v; a = accel(p); v = v + 0.5*dt*a
        if sigma_over_m > 0:
            # SUBCYCLED scattering: nsub sub-steps of dt/nsub per GRAVITY step. This both cuts
            # kappa (and P) proportionally AND gives each sub-step a fresh once-per-step scatter
            # mask, which is what relieves the saturation bias (blocked_frac). Gravity is NOT
            # redone -- only the KDTree is rebuilt per sub-step, which is cheap vs the O(N^2) force.
            if sidm_subcycles == 'auto':
                # measure kappa for a FULL step, then choose nsub so kappa_sub <= KAPPA_TARGET
                from sidm import sidm_kappa
                _k = sidm_kappa(p, v, mass, dt, sigma_over_m, h)
                nsub = int(max(1, math.ceil(_k / KAPPA_TARGET)))
                sidm_diag['kappa_full_step'] = max(sidm_diag.get('kappa_full_step', 0.0), _k)
            else:
                nsub = max(int(sidm_subcycles), 1)
            dts = dt/nsub
            _bf = []
            for _ in range(nsub):
                _d = {}
                v, ns = sidm_scatter(p, v, mass, dts, sigma_over_m, h, rng, diag=_d)
                n_scat += ns
                sidm_diag['n_scatter_total'] = sidm_diag.get('n_scatter_total', 0) + ns
                # per-SUB-STEP safety metrics: these are what must satisfy kappa <= 0.02
                sidm_diag['P_max'] = max(sidm_diag.get('P_max', 0.0), _d.get('P_max', 0.0))
                sidm_diag['kappa'] = max(sidm_diag.get('kappa', 0.0), _d.get('kappa', 0.0))
                _bf.append(_d.get('blocked_frac', 0.0))
                sidm_diag['n_pairs'] = _d.get('n_pairs', 0)
            sidm_diag['blocked_frac'] = max(sidm_diag.get('blocked_frac', 0.0), max(_bf))
            sidm_diag['blocked_frac_mean'] = (
                (sidm_diag.get('blocked_frac_mean', 0.0)*sidm_diag.get('_nbf', 0) + sum(_bf))
                / (sidm_diag.get('_nbf', 0) + len(_bf)))
            sidm_diag['_nbf'] = sidm_diag.get('_nbf', 0) + len(_bf)
            sidm_diag['sidm_subcycles_requested'] = sidm_subcycles
            sidm_diag['sidm_subcycles_used'] = max(sidm_diag.get('sidm_subcycles_used', 0), nsub)
            sidm_diag['subcycles'] = nsub
    return p, v


# ───────────────────────────── EXPERIMENTS ────────────────────────────
# Each returns a dict of headline numbers (and prints a short report).

def exp_observer_gamma(N=200000, seed=1, **kw):
    """Is the apparent inner slope observer-dependent by TRACER ORBIT CLASS? (circularity split)"""
    pos, vel = sample(N, seed, kind='nfw3d')
    eta = circularity(pos, vel)
    r = np.linalg.norm(pos - CEN, axis=1)
    def g_of(mask, lo=0.05, hi=0.40, nb=8):
        rr = r[mask]; e = np.logspace(np.log10(lo), np.log10(hi), nb+1); rc, rho = [], []
        for i in range(nb):
            m = ((rr >= e[i]) & (rr < e[i+1])).sum(); v = 4/3*np.pi*(e[i+1]**3-e[i]**3)
            if m > 30: rc.append((e[i]*e[i+1])**.5); rho.append(m/v)
        return -np.polyfit(np.log(rc), np.log(rho), 1)[0]
    gr, gt = g_of(eta <= 0.4), g_of(eta >= 0.75)
    print(f"  gamma radial={gr:.3f}  tangential={gt:.3f}  Delta={gr-gt:+.3f}"
          f"   {'OBSERVER-DEPENDENT' if abs(gr-gt) >= 0.1 else 'DROP'}")
    return {'gamma_radial': gr, 'gamma_tang': gt, 'dgamma': gr-gt}

def exp_kerr_losscone(N=400000, seed=1, r_g=0.02, **kw):
    """L-null is exactly null for a spherical sink but causal for a Kerr (spin-oriented) loss cone."""
    rng = np.random.default_rng(seed)
    r = 10.0**rng.uniform(-1.0, 0.3, N); f = rng.uniform(0.05, 1.40, N)
    speed = f*np.sqrt(1.0/r)
    ru = rng.normal(size=(N, 3)); ru /= np.linalg.norm(ru, axis=1)[:, None]
    vu = rng.normal(size=(N, 3)); vu /= np.linalg.norm(vu, axis=1)[:, None]
    pos = r[:, None]*ru; vel = speed[:, None]*vu
    E = 0.5*speed**2 - 1.0/r; Lv = np.cross(pos, vel); L = np.linalg.norm(Lv, axis=1)
    ok = (E < -1e-9) & (L > 1e-9); E, L, Lv = E[ok], L[ok], Lv[ok]
    a = -1.0/(2*E); e = np.sqrt(np.clip(1.0 + 2*E*L**2, 0, None)); rp = a*(1-e)
    mu_iso = (Lv/L[:, None]) @ np.array([0.0, 0.0, 1.0])
    def r_mb(chi, mu):
        pro = 2 - chi + 2*np.sqrt(np.clip(1-chi, 0, None)); retro = 2 + chi + 2*np.sqrt(1+chi)
        return retro + (pro - retro)*(mu + 1.0)/2.0
    out = {}
    for chi in (0.0, 0.3, 0.6, 0.9, 0.998):
        cp = float(np.mean(rp < r_g*r_mb(chi, +1.0)/4)); cr = float(np.mean(rp < r_g*r_mb(chi, -1.0)/4))
        out[chi] = cr - cp
        print(f"  chi={chi:5.3f}  cap_pro={cp:.4f} cap_retro={cr:.4f}  asymmetry={cr-cp:+.4f}")
    print(f"  {'BREAKS (L-null causal for Kerr)' if abs(out[0.0])<1e-6 and out[0.998]>0.01 else 'no break'}")
    return {'asym': out}

def exp_stream_shape(qy=0.9, qz=0.7, rc_=0.1, steps=12000, dt=0.02, **kw):
    """Stream-inferred halo shape is a property of the stream's ORBIT FAMILY, not the halo alone."""
    def acc(p):
        D = p[0]**2 + p[1]**2/qy**2 + p[2]**2/qz**2 + rc_**2
        return -np.array([p[0], p[1]/qy**2, p[2]/qz**2])/D
    def run(p0, v0):
        p = np.array(p0, float); v = np.array(v0, float); a = acc(p); pts = np.empty((steps, 3))
        for i in range(steps):
            v = v + 0.5*dt*a; p = p + dt*v; a = acc(p); v = v + 0.5*dt*a; pts[i] = p
        return pts
    fams = {'short-axis tube': ([1.0,0,0.25],[0,1.0,0]), 'long-axis tube': ([0.25,1.0,0],[0,0,1.0]),
            'box / radial': ([0.9,0.5,0.4],[-0.15,0.18,0.15]), 'outer round': ([2.5,0,0],[0,0.85,0.30])}
    bas = {}
    for name, (p0, v0) in fams.items():
        ca, ba = axis_ratios(run(p0, v0)); bas[name] = ba
        print(f"  {name:>18}  b/a(inferred)={ba:.3f}   [c/a={ca:.3f} planarity-confounded]")
    sp = max(bas.values()) - min(bas.values())
    print(f"  truth qy={qy} qz={qz} | inferred in-plane flattening spread={sp:.3f}"
          f"  {'ORBIT-CLASS BIASED' if sp >= 0.05 else 'family-invariant'}")
    return {'b_over_a': bas, 'spread': sp}

def exp_dm_residual(N=200000, seed=1, **kw):
    """Observer-relative DM residual: how much mass a beta-BLIND Jeans observer invents."""
    edges = np.logspace(np.log10(0.05), np.log10(0.80), 22)
    out = {}
    for r_a, lab in [(None,'iso'), (0.40,'0.40'), (0.20,'0.20'), (0.10,'0.10')]:
        pos, vel = sample(N, seed, kind=('isotropic' if r_a is None else 'radial'), r_a=r_a)
        rc, rho, sr2, bet, Mtrue = beta_sigma_profiles(pos, vel, edges)
        M_iso, M_aware = jeans_mass(rc, rho, sr2, bet)
        sel = (rc > 0.10) & (rc < 0.55)
        inv = float(np.median(((M_iso - Mtrue)/Mtrue)[sel])); out[lab] = inv
        print(f"  r_a={lab:>5}  <beta>={np.mean(bet[sel]):+.2f}  "
              f"Jeans-check M_aware/Mtrue={np.median((M_aware/Mtrue)[sel]):.2f}  invented={100*inv:+.0f}%")
    return {'invented': out}

def exp_causal_locus(N=8000, steps=500, theta=0.8, seeds=(1,2,3), **kw):
    """Which orbits -- by CURRENT RADIUS -- causally control the inner slope? (local vs non-local)"""
    shells = [(0.05,0.10),(0.10,0.20),(0.20,0.40),(0.40,0.80)]
    cfg = cfg_for(N); mass = 1.0/N; per = {s: [] for s in shells}; nin = {}
    for sd in seeds:
        pos, vel = sample(N, sd, kind='nfw3d'); r = np.linalg.norm(pos - CEN, axis=1)
        for sh in shells:
            m = (r >= sh[0]) & (r < sh[1]); nin[sh] = int(m.sum())
            g_r = gamma(evolve(pos, radialize(pos, vel, theta, m), cfg, mass, steps)[0])
            g_s = gamma(evolve(pos, sham(pos, vel, theta, m, sd*1000+7), cfg, mass, steps)[0])
            per[sh].append((g_r - g_s)/(nin[sh]/1000.0))
    out = {}
    for sh in shells:
        a = np.array(per[sh]); out[sh] = (a.mean(), a.std(ddof=1))
        print(f"  shell [{sh[0]:.2f},{sh[1]:.2f})  N={nin[sh]:5d}  dgamma/1000 = {a.mean():+.3f} +/- {a.std(ddof=1):.3f}")
    om, os_ = out[shells[-1]]
    print(f"  outer shell significant? {abs(om) > 2*os_}  -> "
          f"{'NON-LOCAL' if abs(om) > 2*os_ else 'LOCAL to r_s'}")
    return {'per_shell': {f'{a}-{b}': out[(a,b)] for a, b in shells}}

def exp_force_law_response(N=8000, steps=500, theta=0.8, seeds=(1,2,3), **kw):
    """Does orbit-reshaping causal RESPONSE discriminate Newton / MOND / superfluid?"""
    cfg = cfg_for(N); mass = 1.0/N; res = {}
    for law in ('newton', 'mond', 'superfluid'):
        rr, tt = [], []
        for sd in seeds:
            pos, vel = sample(N, sd, kind='nfw3d')
            gs = gamma(evolve(pos, sham(pos, vel, theta, None, sd*100+11), cfg, mass, steps, law)[0])
            rr.append(gamma(evolve(pos, radialize(pos, vel, theta), cfg, mass, steps, law)[0]) - gs)
            tt.append(gamma(evolve(pos, tangentialize(pos, vel, theta, None, sd*100+3), cfg, mass, steps, law)[0]) - gs)
        rr, tt = np.array(rr), np.array(tt); A = np.abs(rr) - np.abs(tt)
        res[law] = (A.mean(), A.std(ddof=1))
        print(f"  {law:>11}  radialize={rr.mean():+.3f}+/-{rr.std(ddof=1):.3f}  "
              f"tangential={tt.mean():+.3f}+/-{tt.std(ddof=1):.3f}  asymmetry A={A.mean():+.3f}+/-{A.std(ddof=1):.3f}")
    aN, sN = res['newton']; aS, sS = res['superfluid']
    print(f"  |A(newton)-A(superfluid)| = {abs(aN-aS):.3f} ({abs(aN-aS)/math.sqrt(sN**2+sS**2+1e-9):.1f} sigma)")
    return {'asymmetry': res}

def exp_lyapunov(N=4000, steps=800, delta=1e-5, seeds=(1,), **kw):
    """Cusp irreversibility: exponential (chaotic) vs linear (regular) divergence; lambda(N)."""
    cfg = cfg_for(N); mass = 1.0/N
    pos, vel = sample(N, seeds[0], kind='nfw3d'); vA = radialize(pos, vel, 0.5)
    jr = np.random.default_rng(42)
    accel = make_force('newton', cfg, mass)
    pr, vr = pos.copy(), vA.copy(); ar = accel(pr)
    pp = pos + jr.normal(0, delta, pos.shape); vp = vA.copy(); ap = accel(pp)
    snaps = list(range(0, steps+1, 100)); d = {}
    s = 0
    for tgt in snaps:
        while s < tgt:
            vr = vr + 0.5*DT_DEF*ar; pr = pr + DT_DEF*vr; ar = accel(pr); vr = vr + 0.5*DT_DEF*ar
            vp = vp + 0.5*DT_DEF*ap; pp = pp + DT_DEF*vp; ap = accel(pp); vp = vp + 0.5*DT_DEF*ap
            s += 1
        d[tgt] = float(np.median(np.linalg.norm(pp - pr, axis=1)))
    ts = np.array(snaps[1:])*DT_DEF; amp = np.array([d[t] for t in snaps[1:]])/delta
    lam = np.polyfit(ts, np.log(np.maximum(amp, 1e-30)), 1)[0]
    tcross = RS_DEF/0.8/DT_DEF
    print(f"  N={N}  lambda={lam:.3f}/time   t_lyap={(1/lam)/(tcross*DT_DEF):.1f} t_cross")
    return {'lambda': lam, 'N': N}

def exp_sidm_collapse(N=4000, steps=4000, sigma_over_m=50.0, arm='none', dose=0.0, seeds=(1,),
                      snap=100, rho_r=0.12, sidm_subcycles='auto', **kw):
    """CAMPAIGN: does initial anisotropy causally set SIDM gravothermal collapse time in 3D?

    arm: none | radial | tang | sham | lnull.  Speed-preserving -> positions, per-particle speed,
    KE and virial ratio fixed; only velocity DIRECTIONS change. This separates degree-of-anisotropy
    from the global sigma-profile change that 2510.23705's construction cannot disentangle.
    lnull preserves E, |L| AND pericenter -> a 1D spherical code predicts EXACTLY ZERO response.
    """
    cfg = cfg_for(N); mass = 1.0/N
    out = []
    for sd in seeds:
        if arm == 'om':
            # EQUILIBRIUM Osipkov-Merritt radial IC: beta is SUSTAINED by a stationary DF
            # (unlike a velocity rotation, which makes a non-stationary anisotropy that
            # phase-mixes away in ~1 crossing time). dose = anisotropy radius r_a.
            pos, vel = sample(N, sd, kind='radial', r_a=dose)
            vel0 = vel
        else:
            pos, vel0 = sample(N, sd, kind='isotropic')
        if arm == 'om':
            pass
        elif arm == 'beta':
            # dose is a TARGET BETA; solve the rotation angle for it (map is asymmetric)
            sub_arm, th = dose_for_beta(pos, vel0, dose)
            vel = (vel0 if sub_arm == 'none' else
                   radialize(pos, vel0, th) if sub_arm == 'radial' else
                   tangentialize(pos, vel0, th, seed=sd*7+3))
        elif arm == 'sham':
            # matched null: same rotation MAGNITUDE as the beta arm it controls, random axis
            _, th = dose_for_beta(pos, vel0, dose)
            vel = sham(pos, vel0, th, seed=sd*100+11)
        elif arm == 'lnull':  vel = l_null(pos, vel0, dose)
        else:                 vel = vel0
        b0 = mean_beta(pos, vel)
        rng = np.random.default_rng(sd*13+5); sdiag = {}
        p, v = pos.copy(), vel.copy()
        ts, rcs, cas, bets = [], [], [], []
        for k in range(0, steps, snap):
            p, v = evolve(p, v, cfg, mass, snap, 'newton', sigma_over_m=sigma_over_m, rng=rng,
                          sidm_diag=sdiag, sidm_subcycles=sidm_subcycles)
            r = np.linalg.norm(p - CEN, axis=1)
            ts.append(k+snap)
            rcs.append(float((r < rho_r).sum())/N/(4/3*np.pi*rho_r**3))
            cas.append(axis_ratios(p, rmax=0.5)[0])
            bets.append(mean_beta(p, v))                     # beta(t): isotropization tracking
        ts = np.array(ts); rc_arr = np.array(rcs); bet = np.array(bets)
        imin = int(np.argmin(rc_arr))                        # core forms (rho_c dips) then collapses
        # ROBUST collapse time: FIRST FORWARD crossing of rho_c >= K*rho_c(0).
        # The old "first return above initial AFTER the global minimum" estimator failed two ways:
        #  (a) post-collapse core EVACUATION puts the global min at the END of the run -> NaN for runs
        #      that clearly collapsed (seen at sigma/m=10: rho_c 8.9 -> max 24.3 -> final 1.6);
        #  (b) marginal runs gave spurious late crossings (max rho_c 10.36 vs rho_0 9.62 -> "t=3780").
        # Full rho_c(t) series is always returned so t_collapse can be re-derived at any threshold.
        def _tcol(K):
            j = np.where(rc_arr >= K*rc_arr[0])[0]
            return float(ts[j[0]]) if len(j) else float('nan')
        t_col = _tcol(1.5)
        t_col_k = {f'K{K}': _tcol(K) for K in (1.3, 1.5, 2.0, 3.0)}
        # t_iso = first time |beta| falls below half its initial magnitude (signal-erasure clock)
        t_iso = float('nan')
        if abs(b0) > 0.05:
            j = np.where(np.abs(bet) <= 0.5*abs(b0))[0]
            if len(j): t_iso = float(ts[j[0]])
        out.append(dict(seed=sd, beta_i=b0, beta_f=float(bet[-1]), t_iso=t_iso, t_col_k=t_col_k,
                        sidm_subcycles_requested=sidm_subcycles,
                        sidm_subcycles_used=sdiag.get('sidm_subcycles_used'),
                        kappa_full_step=sdiag.get('kappa_full_step'),
                        P_max=sdiag.get('P_max'), kappa_max=sdiag.get('kappa'),
                        blocked_frac=sdiag.get('blocked_frac'), blocked_frac_mean=sdiag.get('blocked_frac_mean'),
                        n_scatter_total=sdiag.get('n_scatter_total'),
                        rho_min=float(rc_arr[imin]), t_min=float(ts[imin]), t_collapse=t_col,
                        ca_final=float(cas[-1]), beta_series=bet.tolist(), t_series=ts.tolist(),
                        rho_series=rc_arr.tolist()))
        ratio = t_iso/t_col if (t_iso == t_iso and t_col == t_col and t_col > 0) else float('nan')
        print(f"  seed={sd} arm={arm:>6} dose={dose:+.2f} sig={sigma_over_m:g}  "
              f"beta {b0:+.2f}->{bet[-1]:+.2f}  t_iso={t_iso:.0f}  t_collapse={t_col:.0f}  "
              f"t_iso/t_col={ratio:.2f}  c/a={cas[-1]:.3f}")
    return {'arm': arm, 'dose': dose, 'sigma_over_m': sigma_over_m, 'runs': out}


def exp_fingerprints(N=5000, steps=600, seeds=(1,2,3), theta=0.8, **kw):
    """PAPER F3: do orbit / SIDM / feedback occupy SEPARABLE loci in the (df_peri, dgamma) plane?

    Each arm is measured against its OWN matched sham, so the non-equilibrium common mode cancels.
    Ported from the archived collapse_test.py, now with seeds -> error bars (the gap the review flagged).
    """
    cfg = cfg_for(N); mass = 1.0/N; out = {}
    def measure(vel0, pos, evolve_kw, sham_seed):
        pr, vr = evolve(pos, vel0, cfg, mass, steps, **evolve_kw)
        return gamma(pr), f_peri(pr, vr)
    for sd in seeds:
        pos, vel0 = sample(N, sd, kind='nfw3d')
        g_s, f_s = measure(sham(pos, vel0, theta, seed=sd*100+11), pos, dict(law='newton'), sd)
        arms = {}
        for th in (0.4, theta):
            arms[f'radial{th}'] = (radialize(pos, vel0, th), dict(law='newton'))
            arms[f'tang{th}']   = (tangentialize(pos, vel0, th, seed=sd*7+3), dict(law='newton'))
        for som in (10.0, 40.0, 100.0):
            arms[f'sidm{som:.0f}'] = (vel0, dict(law='newton', sigma_over_m=som,
                                                 rng=np.random.default_rng(sd*13+5)))
        for amp in (0.1, 0.2, 0.4):
            arms[f'fb{amp}'] = (vel0, dict(law='newton', feedback_amp=amp))
        for name, (v, ek) in arms.items():
            g, f = measure(v, pos, ek, sd)
            out.setdefault(name, []).append((g - g_s, f - f_s))
    print(f"  {'arm':>10} {'dgamma (mean+/-sd)':>22} {'df_peri (mean+/-sd)':>22}")
    res = {}
    for name, v in out.items():
        a = np.array(v); m, s = a.mean(axis=0), (a.std(axis=0, ddof=1) if len(a) > 1 else np.zeros(2))
        res[name] = dict(dgamma=(float(m[0]), float(s[0])), dfperi=(float(m[1]), float(s[1])))
        print(f"  {name:>10} {m[0]:+10.3f} +/- {s[0]:<7.3f} {m[1]:+10.3f} +/- {s[1]:<7.3f}")
    return res

def exp_scale_flow(N=400000, seed=1, **kw):
    """Does inferred DM flow with the OBSERVER's radial resolution lambda? (ported dm_scaleflow.py)

    Non-circular: sweep the observer's gradient resolution, do NOT smooth the truth (which conserves
    mass and trivially gives zero).
    """
    pos, vel = sample(N, seed, kind='isotropic')
    FE = np.logspace(np.log10(0.04), np.log10(0.9), 60)
    rc, rho, sr2, _, _ = beta_sigma_profiles(pos, vel, FE)
    cf = np.polyfit(np.log(rc), np.log(rho*sr2), 5)          # smooth to isolate curvature from shot noise
    cs = np.polyfit(np.log(rc), np.log(sr2), 4)
    lnP = lambda x: np.polyval(cf, np.log(x))
    s2 = lambda x: np.exp(np.polyval(cs, np.log(x)))
    r = np.linalg.norm(pos - CEN, axis=1); order = np.sort(r)
    Mtrue = lambda x: np.searchsorted(order, x)/N
    def M_inf(rr, lam):
        a, b = rr-lam/2, rr+lam/2
        if a <= rc[0] or b >= rc[-1]: return float('nan')
        return -(rr*s2(rr))*((lnP(b)-lnP(a))/(np.log(b)-np.log(a)))
    LAMS = [0.03, 0.06, 0.12, 0.24]; out = {}
    print(f"  {'r':>6} " + " ".join(f"lam={l:<5}" for l in LAMS) + "   (flow vs finest observer)")
    for rr in (0.15, 0.20, 0.30, 0.40):
        m0 = M_inf(rr, LAMS[0]); row = [M_inf(rr, l)/m0 - 1 if m0 == m0 else float('nan') for l in LAMS]
        out[rr] = row
        print(f"  {rr:6.2f} " + " ".join(f"{100*v:+6.0f}%" if v == v else "   -   " for v in row))
    peak = max(abs(v) for row in out.values() for v in row[1:] if v == v)
    print(f"  peak scale-flow = {100*peak:.0f}%  (vs velocity-access residual of 37-184%)")
    return {'flow': out, 'peak': float(peak)}

def exp_shape_observer(N=20000, steps=1500, r_a=0.10, seed=1, **kw):
    """Is inferred SHAPE observer-dependent by tracer orbit class? (ported shape_observer.py)

    Drives a radial-orbit-unstable halo triaxial, then splits by circularity assigned at t=0.
    """
    cfg = cfg_for(N); mass = 1.0/N
    pos0, vel0 = sample(N, seed, kind='radial', r_a=r_a)
    eta = circularity(pos0, vel0)
    radial, tang = eta <= 0.4, eta >= 0.75
    ca0, _ = axis_ratios(pos0, rmax=0.5)
    pf, _ = evolve(pos0, vel0, cfg, mass, steps, 'newton')
    caf, _ = axis_ratios(pf, rmax=0.5)
    car, _ = axis_ratios(pf[radial], rmax=0.5); cat, _ = axis_ratios(pf[tang], rmax=0.5)
    print(f"  global c/a: {ca0:.3f} -> {caf:.3f}  (triaxiality developed?)")
    print(f"  radial tracers c/a={car:.3f} | tangential c/a={cat:.3f} | Delta={car-cat:+.3f}")
    return {'ca_global': (float(ca0), float(caf)), 'ca_radial': float(car),
            'ca_tang': float(cat), 'delta': float(car-cat)}

EXPERIMENTS = {
    'fingerprints':        exp_fingerprints,
    'scale_flow':          exp_scale_flow,
    'shape_observer':      exp_shape_observer,
    'sidm_collapse':       exp_sidm_collapse,
    'observer_gamma':      exp_observer_gamma,
    'kerr_losscone':       exp_kerr_losscone,
    'stream_shape':        exp_stream_shape,
    'dm_residual':         exp_dm_residual,
    'causal_locus':        exp_causal_locus,
    'force_law_response':  exp_force_law_response,
    'lyapunov':            exp_lyapunov,
}

# recorded published numbers for the regression check (cheap experiments only)
EXPECTED = {
    'kerr_losscone': ('asym', lambda r: abs(r['asym'][0.0]) < 1e-6 and r['asym'][0.998] > 0.15,
                      'chi=0 exactly null; chi->1 asymmetry >0.15'),
    'stream_shape':  ('spread', lambda r: r['spread'] > 0.15, 'in-plane flattening spread >0.15'),
}


# ───────────────────────────────── CLI ────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description='dmlab -- unified DM experiment lab')
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('list')
    r = sub.add_parser('run'); r.add_argument('name'); r.add_argument('--N', type=int)
    r.add_argument('--steps', type=int); r.add_argument('--seeds', type=str)
    r.add_argument('--arm', type=str); r.add_argument('--dose', type=float)
    r.add_argument('--sigma', type=float, dest='sigma_over_m')
    sub.add_parser('verify')
    a = ap.parse_args()

    if a.cmd == 'list':
        for k, f in EXPERIMENTS.items():
            print(f"  {k:<20} {(f.__doc__ or '').strip().splitlines()[0]}")
        return
    if a.cmd == 'verify':
        ok = True
        for name, (key, check, desc) in EXPECTED.items():
            print(f"[{name}] {desc}")
            res = EXPERIMENTS[name]()
            good = check(res); ok &= good
            print(f"  -> {'PASS' if good else 'FAIL'}\n")
        print('REGRESSION:', 'ALL PASS' if ok else 'FAILURES PRESENT')
        sys.exit(0 if ok else 1)
    kw = {}
    if a.N: kw['N'] = a.N
    if a.steps: kw['steps'] = a.steps
    if a.seeds: kw['seeds'] = tuple(int(x) for x in a.seeds.split(','))
    for opt in ('arm', 'dose', 'sigma_over_m'):
        if getattr(a, opt, None) is not None: kw[opt] = getattr(a, opt)
    print(f"[{a.name}] {(EXPERIMENTS[a.name].__doc__ or '').strip()}")
    EXPERIMENTS[a.name](**kw)

if __name__ == '__main__':
    main()
