#!/usr/bin/env python3
"""
nfw_anisotropic_ic.py — near-equilibrium anisotropic NFW initial conditions.
================================================================================
The gating IC for the cusp->core ACCESSIBILITY program (see nbodyv2-real-project
memory): an NFW *cusp* whose velocity anisotropy beta(r) is set by the ANISOTROPIC
(Osipkov-Merritt) JEANS equation, so it sits near virial equilibrium and the
matched-pair f_peri interventions have a stationary starting point.

Why Jeans, not a full DF: a full Eddington/OM distribution function for a SHARPLY
truncated NFW is numerically distorted at the truncation (double-differentiation of
the density blows up at the psi->0 edge -> a sub-virial, contracting halo; measured
Q~0.75-0.92). The anisotropic Jeans equation is a single stable integral, is
guaranteed near-virial (2K=|W| locally), and gives the CORRECT beta(r) by
construction. Its only cost is Gaussian (not exact-DF) velocities -> mild early
relaxation, which is COMMON-MODE across the matched-pair arms (radialize / tangentialize
/ sham-null) and cancels in the differential causal signal. This matches the codebase's
philosophy (nbody_dm_ic.sample_nfw3d already uses isotropic Jeans, Q~1.02).

VALIDATED (2026-07, pass eps = production softening to jeans_sigma_r2 / sample):
  - isotropic: Q~0.99, beta~0, inner-mass HOLDS (+2.4% over 400 steps) -> equilibrium; this is
    the STARTING halo for the cusp->core intervention program.
  - radial: Q~0.98, beta(r) tracks OM 12/12 bins, but inner mass drifts (+22% at r_a=1.5rs).
    Confirmed PHYSICAL radial-orbit instability, not a sampler artifact: drift falls monotonically
    22%->10%->6% as r_a grows 1.5->3->6 rs (weaker anisotropy). The ROI is controllable via r_a
    and is itself part of the science (overshoot<->ROI). For a STABLE radial control use r_a>=6rs.

Anisotropic (OM) Jeans, integrating factor mu(r) = r^2 + r_a^2:
    sigma_r^2(r) = 1/(mu(r) rho(r)) * INT_r^{r_t} mu(r') rho(r') M(<r')/r'^2 dr'
    beta(r)      = r^2 / (r^2 + r_a^2)   (isotropic: r_a -> inf, mu=1, beta=0)
    sigma_t per tangential component = sigma_r * sqrt(1 - beta)
Units G = M_total = 1; x = r/rs; truncated at r_t = c*rs; M(<r) = g(x)/g(c).
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np

Array = np.ndarray


def _g(x):
    return np.log1p(x) - x / (1.0 + x)          # NFW mass function; M(<r) = g(x)/g(c)


class _NFW:
    """Truncated-NFW primitives, code units G = M_total = 1."""
    def __init__(self, rs: float, c: float):
        self.rs, self.c, self.gc, self.rt = rs, c, float(_g(c)), c * rs

    def rho(self, r):
        x = np.asarray(r, float) / self.rs
        rho_s = 1.0 / (4.0 * math.pi * self.rs ** 3 * self.gc)
        with np.errstate(divide="ignore", invalid="ignore"):
            val = rho_s / (x * (1.0 + x) ** 2)
        return np.where((x > 0) & (x < self.c), val, 0.0)

    def Mlt(self, r):
        x = np.clip(np.asarray(r, float) / self.rs, 0.0, self.c)
        return _g(x) / self.gc


def jeans_sigma_r2(nfw: _NFW, r_a: Optional[float], eps: float = 0.0) -> Tuple[Array, Array]:
    """Radial velocity variance sigma_r^2(r) from the anisotropic (OM) Jeans equation.

    eps>0 balances the SOFTENED force the halo is actually evolved under (Plummer
    softening: F(r) = G M(<r) r/(r^2+eps^2)^1.5), so the IC does not re-settle.
    """
    r = np.logspace(math.log10(1e-4 * nfw.rs), math.log10(nfw.rt * (1 - 1e-9)), 8000)
    mu = np.ones_like(r) if r_a is None else (r ** 2 + r_a ** 2)
    force = nfw.Mlt(r) / r ** 2 if eps <= 0 else nfw.Mlt(r) * r / (r ** 2 + eps ** 2) ** 1.5
    integrand = mu * nfw.rho(r) * force                       # G = 1
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)
    tail = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])  # INT_{r_i}^{r_t}
    denom = mu * nfw.rho(r)
    with np.errstate(divide="ignore", invalid="ignore"):
        sig2 = np.where(denom > 0, tail / denom, 0.0)
    return r, np.clip(sig2, 0.0, None)


def sample_nfw_anisotropic(n: int, rs: float, c: float, aniso: str,
                           rng: np.random.Generator, r_a: Optional[float] = None,
                           box_size: Optional[float] = None, eps: float = 0.0) -> Tuple[Array, Array]:
    """(pos, vel), code units G=M=1. aniso in {'isotropic','radial'}; radial needs r_a.

    Pass eps = production softening so the Jeans dispersion balances the softened force.
    """
    if aniso == "isotropic":
        r_a = None
    elif aniso == "radial":
        if r_a is None:
            raise ValueError("radial (Osipkov-Merritt) requires r_a")
    else:
        raise ValueError(aniso)
    nfw = _NFW(rs, c)

    # radii: invert truncated NFW mass CDF  g(x)/g(c) = u
    xg = np.logspace(math.log10(1e-4), math.log10(c), 4096)
    mfrac = _g(xg) / nfw.gc
    u = rng.uniform(mfrac[0], 1.0 - 1e-12, n)
    r = rs * np.interp(u, mfrac, xg)
    rhat = np.empty((n, 3))
    dirs = rng.normal(size=(n, 3))
    rhat = dirs / np.linalg.norm(dirs, axis=1)[:, None]
    pos = r[:, None] * rhat

    # velocities: Gaussian with Jeans sigma_r, sigma_t = sigma_r sqrt(1-beta)
    xg2, sig2g = jeans_sigma_r2(nfw, r_a, eps=eps)
    sig_r = np.sqrt(np.interp(r, xg2, sig2g))
    beta = 0.0 if r_a is None else (r ** 2 / (r ** 2 + r_a ** 2))
    sig_t = sig_r * np.sqrt(np.clip(1.0 - beta, 0.0, None))

    # tangential basis per particle
    ref = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    ref[np.abs(rhat[:, 0]) > 0.9] = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(rhat, ref); e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(rhat, e1)

    vr = rng.normal(0.0, 1.0, n) * sig_r
    vt1 = rng.normal(0.0, 1.0, n) * sig_t
    vt2 = rng.normal(0.0, 1.0, n) * sig_t
    vel = vr[:, None] * rhat + vt1[:, None] * e1 + vt2[:, None] * e2
    vel -= np.mean(vel, axis=0)

    if box_size is not None:
        pos = pos + box_size / 2.0
    return pos, vel


def _beta_profile(pos, vel, center, nbin=12):
    d = pos - center
    r = np.linalg.norm(d, axis=1)
    rhat = d / np.maximum(r[:, None], 1e-30)
    vr = np.sum(vel * rhat, axis=1)
    vt2 = np.sum(vel ** 2, axis=1) - vr ** 2
    edges = np.logspace(np.log10(np.percentile(r, 1)), np.log10(np.percentile(r, 95)), nbin + 1)
    rc, be = [], []
    for i in range(nbin):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if m.sum() < 30:
            continue
        rc.append(math.sqrt(edges[i] * edges[i + 1]))
        be.append(1.0 - np.mean(vt2[m]) / (2.0 * np.mean(vr[m] ** 2)))
    return np.array(rc), np.array(be)


def _minner(pos, center, rthr):
    return float(np.mean(np.linalg.norm(pos - center, axis=1) < rthr))


if __name__ == "__main__":
    from nbody_3d import (SimConfig, integrate, potential_energy_direct, kinetic_energy)
    rs, c = 0.20, 10.0
    cen = np.zeros(3)
    print("NFW anisotropic IC (Jeans) — Q, beta(r), and EVOLVE-AND-HOLD gate\n")

    for aniso, r_a in [("isotropic", None), ("radial", 1.5 * rs)]:
        n = 20000
        rng = np.random.default_rng(7)
        pos, vel = sample_nfw_anisotropic(n, rs, c, aniso, rng, r_a=r_a, eps=0.05)
        sc = SimConfig(model="direct_isolated", integrator="leapfrog_kdk", init="nfw3d",
                       seed=1, n=n, eps=0.05, box_size=2.0, plummer_a=rs)
        mass = 1.0 / n
        Q0 = 2.0 * kinetic_energy(vel, mass) / abs(potential_energy_direct(pos, mass, sc))
        rc, be = _beta_profile(pos, vel, cen)
        label = aniso + (f" (r_a={r_a:.2f})" if r_a else "")
        print(f"[{label}]  Q0 = {Q0:.3f}")
        bad = 0
        for r0, b0 in zip(rc, be):
            exp = (r0 ** 2 / (r0 ** 2 + r_a ** 2)) if r_a else 0.0
            if abs(b0 - exp) > 0.12 and r0 > 0.06:
                bad += 1
        print(f"   beta(r) matches OM within 0.12 in {len(rc)-bad}/{len(rc)} bins (r>0.06)")

        # EVOLVE-AND-HOLD gate: does the inner mass fraction stay put?
        ne = 4000
        rng2 = np.random.default_rng(11)
        pe, ve = sample_nfw_anisotropic(ne, rs, c, aniso, rng2, r_a=r_a, eps=0.05)
        sc2 = SimConfig(model="direct_isolated", integrator="leapfrog_kdk", init="nfw3d",
                        seed=1, n=ne, eps=0.05, box_size=2.0, plummer_a=rs, steps=400, dt=0.005)
        m0 = _minner(pe, cen, 0.2)
        snaps = integrate(pe, ve, 1.0 / ne, sc2, [200, 400], use_numba=False)
        mf = _minner(snaps[400][0], cen, 0.2)
        drift = (mf - m0) / m0
        print(f"   M(<0.2) fraction: t0={m0:.3f} -> t400={mf:.3f}  (drift {drift*100:+.1f}%)  "
              f"{'HOLDS' if abs(drift) < 0.05 else 'DRIFTS -> not equilibrium'}\n")
    print("GATE: Q~1, beta tracks OM, and |inner-mass drift| < 5% over 400 steps.")
