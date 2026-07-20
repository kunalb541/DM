#!/usr/bin/env python3
"""
nbody_dm_ic.py — dark-matter (NFW) initial conditions for the N-body battery.
================================================================================
Tier-A of the DM upgrade (see dark-matter-reframe memory / PREREG_collisionality_
crossover.md): drop a *real* dark-matter halo profile into the existing battery so
both the predictive (nbody_stress.run_stress) and causal (nbody_aws_battery.
_cell_worker) tracks can run on NFW haloes with NO change to the observable or
intervention machinery.

Profiles
--------
  nfw3d     — cuspy NFW halo (rho ~ r^-1 inner, the canonical CDM cusp).
              Radii sampled from the exact truncated NFW enclosed-mass CDF;
              velocities from the isotropic Jeans equation (near-equilibrium,
              validated by the virial ratio Q = 2K/|W| ~ 1 self-test below).

Cusp <-> core framing: NFW (nfw3d, cuspy) vs Plummer (plummer3d, cored) are the
two ends of the cusp-core dichotomy the causal f_peri handle is meant to control.
Both plug into the same StressConfig machinery. A cored-NFW / Burkert variant is
a natural next addition (see _EXTENSION note at bottom).

Units: code units G = M_total = 1, matching the rest of the codebase. Scale
radius rs is taken from cfg.plummer_a (reused as the halo scale length); the
concentration c is cfg.nfw_c (added to StressConfig with default 10.0). The halo
is truncated at the virial radius r_vir = c * rs and centred at box_size/2.
"""
from __future__ import annotations
import math
from typing import Tuple

import numpy as np

Array = np.ndarray

# Grid resolution for the NFW inverse-CDF and Jeans integral. 4096 log-spaced
# nodes gives sub-0.1% interpolation error on both; cheap (done once per sample).
_NGRID = 4096
_XMIN = 1e-4   # inner x=r/rs cutoff for the velocity integral (avoids log divergence)


def _sphere_directions(rng: np.random.Generator, n: int) -> Array:
    """n unit vectors uniformly on S^2 (matches nbody_3d._sphere_directions)."""
    cth = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * math.pi, n)
    sth = np.sqrt(1.0 - cth * cth)
    return np.column_stack([sth * np.cos(phi), sth * np.sin(phi), cth])


def _g(x: Array | float) -> Array | float:
    """NFW dimensionless mass function g(x) = ln(1+x) - x/(1+x).  M(<r) = g(x)/g(c)."""
    return np.log1p(x) - x / (1.0 + x)


def _rho_tilde(x: Array) -> Array:
    """NFW density shape (without rho_s): 1 / [ x (1+x)^2 ]."""
    return 1.0 / (x * (1.0 + x) ** 2)


def nfw_concentration(cfg) -> float:
    return float(getattr(cfg, "nfw_c", 10.0))


def _jeans_sigma2_of_x(c: float, rs: float) -> Tuple[Array, Array]:
    """Isotropic-Jeans 1-D velocity variance sigma^2(x) on a log x-grid in (0, c].

    Isotropic Jeans (no anisotropy): d(rho sigma_r^2)/dr = -rho G M(<r)/r^2, so
        sigma_r^2(r) = (1/rho(r)) * INT_r^{r_t} rho(r') G M(<r')/r'^2 dr'.
    In dimensionless x = r/rs with G = M_tot = 1, M(<r') = g(x')/g(c):
        sigma^2(x) = 1/(rs g(c) rho~(x)) * INT_x^{c} rho~(x') g(x') / x'^2 dx'.
    Truncated at the virial radius r_t = c rs (isolated halo).  Returns
    (x_grid, sigma2_grid) with sigma2 the 1-D (per-component) variance.
    """
    x = np.logspace(math.log10(_XMIN), math.log10(c), _NGRID)
    integrand = _rho_tilde(x) * _g(x) / x ** 2
    # Cumulative integral from the OUTER edge inward: INT_x^c ... dx'
    # trapezoid on the (increasing) x-grid, then reverse-accumulate.
    dseg = np.diff(x)
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * dseg          # per-interval
    tail = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])   # INT_{x_i}^{c}
    sigma2 = tail / (rs * _g(c) * _rho_tilde(x))
    sigma2 = np.clip(sigma2, 0.0, None)
    return x, sigma2


def sample_nfw3d(rng: np.random.Generator, cfg) -> Tuple[Array, Array]:
    """Truncated NFW halo: exact enclosed-mass radii + isotropic-Jeans velocities.

    Returns (pos, vel) in code units, centred at box_size/2, total mass 1
    (per-particle mass 1/N is applied by the integrator, not here).
    """
    n = int(cfg.n)
    rs = float(cfg.plummer_a)
    c = nfw_concentration(cfg)
    box = float(cfg.box_size)

    # ── radii: invert the truncated NFW mass CDF  g(x)/g(c) = u  ────────────────
    x_grid = np.logspace(math.log10(_XMIN), math.log10(c), _NGRID)
    mfrac = _g(x_grid) / _g(c)                     # monotone 0->1
    u = rng.uniform(mfrac[0], 1.0 - 1e-12, n)      # avoid the x->0 singular tail
    x = np.interp(u, mfrac, x_grid)
    r = rs * x
    center = np.full(3, box / 2.0)
    pos = center + r[:, None] * _sphere_directions(rng, n)

    # ── velocities: isotropic Jeans sigma(r), Gaussian per component  ───────────
    xg, sig2g = _jeans_sigma2_of_x(c, rs)
    sig2 = np.interp(x, xg, sig2g)
    sigma = np.sqrt(np.clip(sig2, 0.0, None))
    vel = rng.normal(0.0, 1.0, (n, 3)) * sigma[:, None]
    vel -= np.mean(vel, axis=0)
    return pos, vel


# ── self-test: sample an NFW halo and report the initial virial ratio ──────────
def _selftest() -> None:
    """python nbody_dm_ic.py — sample NFW at a few N and print Q = 2K/|W|.

    A near-equilibrium isotropic-Jeans halo should give Q ~ 1. Large departures
    would flag a bug in the sampler or Jeans integral before any battery run.
    """
    from dataclasses import dataclass

    @dataclass
    class _Cfg:
        n: int
        plummer_a: float = 0.20
        nfw_c: float = 10.0
        box_size: float = 2.0
        G: float = 1.0
        eps: float = 0.05

    try:
        from nbody_3d import SimConfig, potential_energy_direct, kinetic_energy
        have_pe = True
    except Exception:
        have_pe = False

    print(f"{'N':>7}{'c':>6}{'rs':>7}{'r95/box':>10}{'Q=2K/|W|':>12}")
    for n in (1024, 2048, 4096):
        cfg = _Cfg(n=n)
        rng = np.random.default_rng(12345)
        pos, vel = sample_nfw3d(rng, cfg)
        r = np.linalg.norm(pos - cfg.box_size / 2.0, axis=1)
        r95 = float(np.percentile(r, 95.0)) / cfg.box_size
        if have_pe:
            sc = SimConfig(model="direct_isolated", integrator="leapfrog_kdk",
                           init="nfw3d", seed=1, n=n, eps=cfg.eps,
                           box_size=cfg.box_size, plummer_a=cfg.plummer_a)
            mass = 1.0 / n
            ke = kinetic_energy(vel, mass)
            pe = potential_energy_direct(pos, mass, sc)
            Q = 2.0 * ke / abs(pe) if pe != 0 else float("nan")
        else:
            Q = float("nan")
        print(f"{n:>7}{cfg.nfw_c:>6.1f}{cfg.plummer_a:>7.2f}{r95:>10.3f}{Q:>12.3f}")
    print("\nQ ~ 1 => near virial equilibrium (isotropic-Jeans NFW is not an exact "
          "DF, so a mild offset and early relaxation are expected, as with the "
          "codebase's Hernquist/Plummer samplers).")


# _EXTENSION: cored-NFW / Burkert — same pattern. Replace _rho_tilde and _g with
# the cored profile's density + mass function; the sampler and Jeans integral are
# otherwise identical. Gives a matched cuspy/cored pair for the cusp-core causal test.

if __name__ == "__main__":
    _selftest()
