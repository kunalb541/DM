#!/usr/bin/env python3
"""
feedback.py — baryonic-feedback core formation via repeated non-adiabatic potential fluctuation
(the Pontzen & Governato 2012 mechanism, in controlled toy form).

A central "gas" mass oscillates in amplitude (gas repeatedly driven out by supernovae and
re-accreted). When the oscillation period is <~ the local orbital time, dark-matter orbits
cannot follow adiabatically and gain energy each cycle -> the cusp heats and flattens. This
is GENTLE and physical, unlike the earlier velocity-kick toy (which just ejected particles;
that version is retained as feedback_kick for reference but is NOT used).

DM acceleration each step = gravity(DM) + gas potential accel, with
  M_gas(t) = M_amp * (0.5 - 0.5 cos(2 pi t / T_osc))            (oscillates 0 -> M_amp -> 0)
  a_gas_i  = -G M_gas(t) (x_i - c) / (|x_i-c|^2 + r_gas^2)^{3/2}   (Plummer-softened central mass)
Dose axes: M_amp (gas mass, = strength) and T_osc (period, = frequency). Energy is injected by
design. VALIDATION-FIRST, resolved metrics (r>eps). Units G = M_total = 1.
"""
from __future__ import annotations
import types, numpy as np
from nbody_3d import acceleration, apply_pbc


def feedback_kick(pos, vel, center, r_burst, v_kick, rng):
    """DEPRECATED blunt model (ejects particles); kept for reference only."""
    vel = vel.copy(); r = np.linalg.norm(pos - center, axis=1); inner = r < r_burst
    n = int(inner.sum())
    if n and v_kick > 0:
        vel[inner] += rng.normal(0.0, v_kick, (n, 3))
    return vel, n


def _gas_acc(pos, center, M_gas, r_gas):
    dx = pos - center
    r2 = np.sum(dx * dx, axis=1) + r_gas ** 2
    return -M_gas * dx / (r2[:, None] ** 1.5)


def evolve_feedback(pos, vel, mass, cfg, steps, M_amp, r_gas, T_osc, snap_every=50):
    """Leapfrog + oscillating central gas potential (Pontzen-Governato mechanism)."""
    center = np.full(3, cfg.box_size / 2.0)
    pos, vel = pos.copy(), vel.copy()
    dt = cfg.dt
    def total_acc(p, step):
        Mg = M_amp * (0.5 - 0.5 * np.cos(2 * np.pi * step / T_osc)) if (M_amp > 0 and T_osc > 0) else 0.0
        a = acceleration(p, mass, cfg, use_numba=True)
        return a + (_gas_acc(p, center, Mg, r_gas) if Mg != 0.0 else 0.0)
    acc = total_acc(pos, 0)
    snaps = {0: (pos.copy(), vel.copy())}
    for step in range(1, steps + 1):
        vel = vel + 0.5 * dt * acc
        pos = apply_pbc(pos + dt * vel, cfg)
        acc = total_acc(pos, step)
        vel = vel + 0.5 * dt * acc
        if step % snap_every == 0:
            snaps[step] = (pos.copy(), vel.copy())
    return snaps


if __name__ == "__main__":
    from nbody_dm_ic import sample_nfw3d
    from nbody_3d import SimConfig
    N, RS, EPS, STEPS = 3000, 0.20, 0.05, 800
    R_GAS, T_OSC, SEEDS = 0.10, 100, [1, 2, 3]     # T_osc=100 steps ~ half a central dyn.time (non-adiabatic)
    AMPS = [0.0, 0.05, 0.10, 0.20]                  # gas mass as fraction of total (M(<0.2)~0.13)
    cen = np.full(3, 1.0)
    icfg = types.SimpleNamespace(n=N, plummer_a=RS, nfw_c=10.0, box_size=2.0)
    sc = SimConfig(model='direct_isolated', integrator='leapfrog_kdk', init='nfw3d',
                   seed=1, n=N, eps=EPS, box_size=2.0, plummer_a=RS, steps=STEPS, dt=0.005)

    def mf(p, r): return float(np.mean(np.linalg.norm(p - cen, axis=1) < r))
    def slope(p):
        r = np.linalg.norm(p - cen, axis=1); e = np.array([0.06, 0.10, 0.16, 0.25]); rc, rho = [], []
        for i in range(len(e) - 1):
            m = ((r >= e[i]) & (r < e[i + 1])).sum(); v = 4 / 3 * np.pi * (e[i + 1] ** 3 - e[i] ** 3)
            if m > 5:
                rc.append((e[i] * e[i + 1]) ** .5); rho.append(m / v)
        return -np.polyfit(np.log(rc), np.log(rho), 1)[0] if len(rc) > 2 else float('nan')

    print(f"FEEDBACK (potential-fluctuation) DOSE-RESPONSE  (N={N}, {STEPS} steps, T_osc={T_OSC}, r_gas={R_GAS})")
    print(f"{'M_gas/Mtot':>11} {'inner slope':>16} {'M(<0.1) [x1e-3]':>18}")
    res = {}
    for amp in AMPS:
        sl, mm = [], []
        for s in SEEDS:
            pos, vel = sample_nfw3d(np.random.default_rng(s), icfg)
            snaps = evolve_feedback(pos, vel, 1.0 / N, sc, STEPS, amp, R_GAS, T_OSC)
            pf = snaps[STEPS][0]; sl.append(slope(pf)); mm.append(mf(pf, 0.1) * 1e3)
        sl, mm = np.array(sl), np.array(mm)
        res[amp] = (np.nanmean(sl), np.nanstd(sl) / np.sqrt(len(sl)))
        print(f"{amp:11.2f}   {np.nanmean(sl):.2f} +- {np.nanstd(sl)/np.sqrt(len(sl)):.2f}      {mm.mean():6.1f} +- {mm.std()/np.sqrt(len(mm)):.1f}")
    print("\nPASS if amp=0 keeps cusp and amp>0 flattens it GENTLY (a core, not an evacuated hole).")
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        SP = "/private/tmp/claude-501/-Users-kunalbhatia-Desktop-Research-nbody/aa638959-d922-48ba-80e3-6df794ae9f7a/scratchpad"
        fig, ax = plt.subplots(figsize=(6.0, 4.0))
        ax.errorbar(AMPS, [res[a][0] for a in AMPS], yerr=[res[a][1] for a in AMPS], fmt='o-', capsize=3, color='C4')
        ax.axhline(res[0.0][0], ls='--', c='grey', lw=.8, label='CDM (no feedback)')
        ax.set_xlabel('gas mass amplitude  M_gas / M_tot'); ax.set_ylabel('inner density slope (cusp~1.5, core~0)')
        ax.set_title('Feedback (potential fluctuation): cusp flattening vs gas mass')
        ax.legend(fontsize=9); fig.tight_layout(); fig.savefig(f'{SP}/fig_feedback_dose.png', dpi=140)
        print('fig_feedback_dose saved.')
    except Exception as e:
        print('plot skipped:', e)
