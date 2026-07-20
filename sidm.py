#!/usr/bin/env python3
"""
sidm.py — self-interacting dark matter via probabilistic pairwise Monte-Carlo scattering.

The #1 rival to orbit-reshaping for making cores. Adds a small chance per step for nearby
particles to elastically scatter (isotropic, momentum- and energy-conserving), on top of
gravity. VALIDATION-FIRST: this module is only trustworthy once it (a) forms a core on an
isolated NFW at sigma/m>0, and (b) recovers CDM (no core) at sigma/m=0. The __main__ block
runs exactly that gate. Do NOT use in the mechanism scoreboard until the gate passes.

Method (simplified Rocha 2013 / Robertson 2017): for pairs within smoothing length h,
P_scatter = (sigma/m) * rho_local * |v_rel| * dt; on scatter, redraw the relative-velocity
direction isotropically about the pair centre of mass (|v_rel| preserved => elastic).
Units: G = M_total = 1 (code units); sigma/m is in code units (physical cm^2/g is a
calibration set by the halo's physical mass+size).
"""
from __future__ import annotations
import numpy as np


def sidm_scatter(pos, vel, mass, dt, sigma_over_m, h, rng, diag=None):
    """One elastic-scattering sub-step (in place on a vel copy). Returns (vel, n_scatter).

    If `diag` is a dict, records SAFETY diagnostics required to certify the Monte-Carlo:
      P_max        largest single-pair scatter probability this step (must be << 1)
      kappa        expected scatters per particle this step = dt/t_scat (K&S: <=0.02 for O(10%))
      n_scatter    actual scatters executed
      n_pairs      candidate pairs in the kernel
      blocked_frac fraction of candidate pairs SKIPPED because a partner had already scattered
                   -- this is the once-per-step saturation bias, which grows as kappa rises
    """
    from scipy.spatial import cKDTree
    vel = vel.copy()
    if sigma_over_m <= 0:
        return vel, 0
    tree = cKDTree(pos)
    pairs = tree.query_pairs(h, output_type='ndarray')
    if len(pairs) == 0:
        return vel, 0
    Vh = (4.0 / 3.0) * np.pi * h ** 3          # kernel volume for the local density estimate
    rho = mass / Vh                             # density contributed by one partner in the kernel
    scattered = np.zeros(len(pos), bool)
    n_sc = 0
    _pmax = 0.0; _psum = 0.0; _blocked = 0
    for idx in rng.permutation(len(pairs)):
        i, j = pairs[idx]
        if scattered[i] or scattered[j]:
            _blocked += 1
            continue
        dv = vel[i] - vel[j]
        vrel = np.sqrt(dv @ dv)
        P = sigma_over_m * rho * vrel * dt
        if P > _pmax: _pmax = P
        _psum += P
        if rng.random() < P:
            vcm = 0.5 * (vel[i] + vel[j])
            n = rng.normal(size=3); n /= np.sqrt(n @ n)
            vel[i] = vcm + 0.5 * vrel * n
            vel[j] = vcm - 0.5 * vrel * n
            scattered[i] = scattered[j] = True
            n_sc += 1
    if diag is not None:
        diag['P_max'] = float(_pmax)
        diag['kappa'] = float(2.0*_psum/len(pos))      # expected scatters per particle per step
        diag['n_scatter'] = int(n_sc)
        diag['n_pairs'] = int(len(pairs))
        diag['blocked_frac'] = float(_blocked/max(len(pairs), 1))
    return vel, n_sc


def evolve_sidm(pos, vel, mass, cfg, steps, sigma_over_m, h, rng, snap_every=50):
    """Leapfrog KDK gravity + SIDM scatter (operator-split), returning {step:(pos,vel)}."""
    from nbody_3d import acceleration, apply_pbc
    pos, vel = pos.copy(), vel.copy()
    acc = acceleration(pos, mass, cfg, use_numba=True)
    dt = cfg.dt
    snaps = {0: (pos.copy(), vel.copy())}
    tot_sc = 0
    for step in range(1, steps + 1):
        vel = vel + 0.5 * dt * acc
        pos = apply_pbc(pos + dt * vel, cfg)
        acc = acceleration(pos, mass, cfg, use_numba=True)
        vel = vel + 0.5 * dt * acc
        vel, nsc = sidm_scatter(pos, vel, mass, dt, sigma_over_m, h, rng)  # scatter sub-step
        tot_sc += nsc
        if step % snap_every == 0:
            snaps[step] = (pos.copy(), vel.copy())
    return snaps, tot_sc


if __name__ == "__main__":
    import types, math
    from nbody_dm_ic import sample_nfw3d
    from nbody_3d import SimConfig, kinetic_energy, potential_energy_direct
    N, RS, EPS = 3000, 0.20, 0.05
    STEPS, H = 600, 0.10
    cen = np.full(3, 1.0)

    def mfrac(p, r):
        return float(np.mean(np.linalg.norm(p - cen, axis=1) < r))

    print(f"SIDM VALIDATION GATE  (isolated NFW, N={N}, {STEPS} steps, h={H})")
    print("Does a core form for sigma/m>0 and NOT for sigma/m=0?\n")
    print(f"{'sigma/m':>8} {'scatters':>10} {'M(<0.05) t0->tf':>20} {'change':>9} {'energy drift':>13}")
    icfg = types.SimpleNamespace(n=N, plummer_a=RS, nfw_c=10.0, box_size=2.0)
    sc = SimConfig(model="direct_isolated", integrator="leapfrog_kdk", init="nfw3d",
                   seed=1, n=N, eps=EPS, box_size=2.0, plummer_a=RS, steps=STEPS, dt=0.005)
    mass = 1.0 / N
    for som in [0.0, 10.0, 40.0]:
        rng = np.random.default_rng(1)
        pos, vel = sample_nfw3d(np.random.default_rng(1), icfg)
        E0 = kinetic_energy(vel, mass) + potential_energy_direct(pos, mass, sc)
        snaps, nsc = evolve_sidm(pos, vel, mass, sc, STEPS, som, H, rng)
        pf, vf = snaps[STEPS]
        Ef = kinetic_energy(vf, mass) + potential_energy_direct(pf, mass, sc)
        m0, mf = mfrac(pos, 0.05), mfrac(pf, 0.05)
        print(f"{som:8.1f} {nsc:10d}   {m0:.3f} -> {mf:.3f}     {100*(mf-m0)/m0:+7.1f}%  {abs((Ef-E0)/E0):13.2e}")
    print("\nPASS if: sigma/m=0 keeps M(<0.05) ~flat (CDM), and sigma/m>0 REDUCES it (core forms),")
    print("more so at larger sigma/m, with energy drift small. Then SIDM is trustworthy for the scoreboard.")
