"""Tests for the live lab. Run: python -m pytest tests/ -q

These assert the INVARIANTS the causal-intervention method depends on. If any of these break,
every published differential in this repo is invalid.
"""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dmlab as D

N = 2000

def _halo():
    return D.sample(N, 1, kind='nfw3d')

# --- the interventions must be speed-preserving (this is what makes them matched pairs) ---
def test_rotation_is_exactly_speed_preserving_before_momentum_rezero():
    """The ROTATION itself is exact. _apply() then subtracts the mean velocity to re-zero
    momentum, and THAT is what perturbs individual speeds (|<v>| ~ 0.008 vs speed ~ 0.66).
    Measured: rotation alone gives median |dv|/v = 0 exactly."""
    pos, vel = _halo()
    rhat, v_r, vt, v_t, sp = D._frame(pos, vel)
    that = vt / np.maximum(v_t, 1e-12)[:, None]
    phi = np.arctan2(v_t, v_r)
    phi2 = np.where(phi <= math.pi/2, np.maximum(phi-0.6, 0.0), np.minimum(phi+0.6, math.pi))
    raw = sp[:, None]*(np.cos(phi2)[:, None]*rhat + np.sin(phi2)[:, None]*that)
    assert np.allclose(np.linalg.norm(raw, axis=1), sp, atol=1e-12)

def test_KE_preserved_to_1e3_after_momentum_rezero():
    """KE is held to ~1 part in 1e4 (measured), NOT exactly -- the momentum re-zero costs that.
    It is common-mode across real and sham arms, so matched-pair differentials are unaffected.
    The paper must say 'to O(1e-4)', not 'exactly by construction'."""
    pos, vel = _halo()
    for f in (D.radialize, D.tangentialize, D.sham):
        assert abs(np.sum(f(pos, vel, 0.6)**2)/np.sum(vel**2) - 1) < 1e-3

def test_sham_beta_shift_within_estimator_noise():
    """Sham must be anisotropy-neutral WITHIN the noise of the beta estimator itself, which is
    ~0.05-0.1 at these N (an 'isotropic' IC reads beta=+0.19 at N=2000 from sampling noise alone).
    So compare the shift against the seed-to-seed scatter of beta, not against zero."""
    shifts, base = [], []
    for sd in (1, 2, 3, 4, 5):
        pos, vel = D.sample(4000, sd, kind='nfw3d')
        b0 = D.mean_beta(pos, vel); base.append(b0)
        shifts.append(D.mean_beta(pos, D.sham(pos, vel, 0.8)) - b0)
    noise = np.std(base, ddof=1)
    assert abs(np.mean(shifts)) < 3*max(noise, 0.03), \
        f"sham shift {np.mean(shifts):+.3f} exceeds 3x beta-estimator noise {noise:.3f}"

def test_l_null_preserves_L_and_pericenter():
    """l_null is the instrument that is null for spherical systems: it must hold |L| AND r_peri."""
    pos, vel = _halo()
    v2 = D.l_null(pos, vel, 0.8)
    L1 = np.linalg.norm(np.cross(pos - D.CEN, vel), axis=1)
    L2 = np.linalg.norm(np.cross(pos - D.CEN, v2), axis=1)
    assert abs(np.median(L2/np.maximum(L1, 1e-12)) - 1) < 1e-3
    rp1, rp2 = D.pericenters(pos, vel), D.pericenters(pos, v2)
    assert abs(np.median(rp2/np.maximum(rp1, 1e-12)) - 1) < 1e-2

def test_momentum_rezeroed():
    pos, vel = _halo()
    for f in (D.radialize, D.tangentialize, D.sham, D.l_null):
        assert np.linalg.norm(np.mean(f(pos, vel, 0.5), axis=0)) < 1e-9

# --- measurement sanity ---
def test_gamma_recovers_nfw_cusp():
    pos, _ = _halo()
    g = D.gamma(pos)
    assert 0.8 < g < 2.5, f"NFW inner slope out of range: {g}"

def test_dose_for_beta_hits_target():
    """The dose->beta map is strongly asymmetric; the solver must hit the requested beta."""
    pos, vel = D.sample(3000, 1, kind='isotropic')
    for target in (0.3, -0.3):
        arm, th = D.dose_for_beta(pos, vel, target)
        v2 = D.radialize(pos, vel, th) if arm == 'radial' else D.tangentialize(pos, vel, th)
        assert abs(D.mean_beta(pos, v2) - target) < 0.08

def test_om_ic_is_radially_anisotropic():
    pos, vel = D.sample(3000, 1, kind='radial', r_a=0.18)
    assert D.mean_beta(pos, vel) > 0.25

# --- published-result regressions (same checks as `dmlab.py verify`) ---
def test_kerr_losscone_null_then_causal():
    r = D.exp_kerr_losscone(N=100000)
    assert abs(r['asym'][0.0]) < 1e-6, "spherical sink must be EXACTLY null"
    assert r['asym'][0.998] > 0.10, "Kerr sink must be causal"

def test_stream_shape_is_orbit_class_biased():
    r = D.exp_stream_shape(steps=4000)
    assert r['spread'] > 0.05
