> **STATUS: CANONICAL PLAN.**
> Near-term paper = the SIDM anisotropy result below. `PAPER_PLAN.md` (causal fingerprints)
> is FUTURE WORK until its F3 battery has seeded error bars.

## Near-term paper (canonical)
**"A matched-intervention test of velocity anisotropy in SIDM gravothermal collapse"**

**Thesis.** Kamionkowski & Sigurdson (2025) report >2x variation in SIDM gravothermal collapse
time with anisotropic velocity distributions, from a 1D spherically-symmetric solver. We apply a
3D matched-pair causal intervention that changes velocity **direction** while holding the density
field, positions, and per-particle speeds fixed (KE to O(1e-4)). This isolates *anisotropy* from
the distribution-function and sigma-profile changes that co-vary in self-consistent anisotropic ICs.

**Scope — deliberately narrow.** We test whether **INITIAL** anisotropy causally changes collapse
time. We do **not** claim to hold anisotropy fixed during collapse: we measure that it is erased
fast (~1 crossing time by phase mixing for a rotated IC; ~20 steps by scattering for an equilibrium
OM IC), and t_iso/t_collapse ~ 0.005-0.019 across sigma/m = 5-80. **That erasure timescale is a
RESULT, reported prominently — not a nuisance.** It is also common to the 1D work (their
t_scat ~ 110 Myr vs ~16 Gyr collapse), so any anisotropy effect must be an EARLY IMPRINT on the
density/energy structure rather than a sustained state.

**Novelty claim.** NOT "anisotropy matters" (known). The claim is: *a matched causal intervention
isolates the anisotropy/orbit-channel effect from distribution-function and sigma-profile confounds*,
and bounds it in 3D.

**Target:** MNRAS main journal.

**BLOCKING GATE:** no production compute until the SIDM dt-convergence test passes (see
`fleet/conv_cell.py`). If t_collapse or the rho trajectory shifts >10-15% between dt=0.005 and
dt=0.00125, all existing SIDM claims become calibration-only and must be rerun with smaller dt or
subcycling before Phase 1.

---

# Paper plan + full-scale AWS plan (for external review)

Every number below is traceable to a measurement in this repo. Provenance is given as
`[run-id / script]`. Claims I previously made and later **refuted** are listed in §5 so they are
not re-proposed. Written after 103 AWS cells across cal1–cal6 (~$4 total, 0 orphaned instances).

---

## 1. VERIFIED RESULTS (survive the fixed collapse estimator, threshold-robust)

| # | result | numbers | provenance |
|---|--------|---------|------------|
| V1 | **SIDM collapse is GRAVOTHERMAL, not two-body relaxation** | d(ln t_col)/d(ln N) = +0.01 / −0.08 / −0.06 at K=1.3/1.5/2.0, over N=2000→8000, while t_relax grows 3.4× | cal6 `F_N*` |
| V2 | **Matched-pair sham control is sound (scatter, not bias)** | none vs sham: 1.6σ / 0.2σ / 0.5σ at K=1.3/1.5/2.0 | cal6 `F_none/F_sham` |
| V3 | **SIDM erases velocity anisotropy 1–2 orders of magnitude faster than it collapses** | t_iso/t_collapse ≈ 0.005–0.019 across σ/m = 5→80 | cal6 `F_sig*` |
| V4 | **The erasure channel is COLLISIONAL, not collisionless** | equilibrium OM β retained 61% at t=1000 under pure gravity (σ/m=0) vs **5%** under SIDM; rotated β retains 34% vs 4% | cal4 `D_om_grav/D_om_sidm` |
| V5 | **A factor-2 β-dependence of collapse time is EXCLUDED in 3D** | isotropic 1860±20 vs radial groups 1985–2707; a 2× effect would require ~3720 | cal5 `E_*` |
| V6 | **Run-to-run scatter is DYNAMICAL, not Poisson** | sd/mean = 16.2 / 10.2 / 9.4 / 13.5 % at N = 2k/4k/6k/8k → N^−0.21 (Poisson −0.5) | cal6 `F_N*` |

**Bounded null (V5 detail):** a ≲30% β-effect is **NOT** excluded at current power. The apparent
trend's SIGN FLIPS with the collapse threshold (slope +396 / +711 / −919 at K=1.3/1.5/2.0) — the
signature of noise, not signal.

---

## 2. PAPER 1 — ready now, does not depend on any further compute
**"Observer-dependent dark matter: the inferred halo is a property of the tracer, not the halo"**

Five confirmed, independently-measured effects on one axis (velocity/orbit structure), each
reframing a known degeneracy as a *measurable observer-relative residual*:

| effect | magnitude | script |
|---|---|---|
| inner slope γ depends on tracer orbit class (circularity split) | **Δγ = 0.24** | `dmlab observer_gamma` |
| inferred triaxiality depends on tracer orbit class | **Δ(c/a) = 0.44** | `dmlab run shape_observer` |
| stream-inferred shape is orbit-family biased | in-plane b/a **0.71→0.90** for the same halo | `dmlab stream_shape` |
| β-blind Jeans observer invents mass | **+37% → +184%** as β 0.28→0.80; isotropic control ≈0 | `dmlab dm_residual` |
| observer-resolution scale-flow | ~10% (sub-dominant axis) | `dmlab run scale_flow` |

**Honest framing (this is what fixes the two prior rejections):** we are NOT claiming to discover
the mass–anisotropy degeneracy (Binney; Wolf+2010) or finite-resolution Jeans bias. The
contribution is (a) defining an observer-relative residual fraction and (b) measuring it across
density, geometry, and streams with exact ground truth, showing it is **dominantly a
velocity-information phenomenon (37–184%)** and only weakly a spatial-resolution one (~10%).
Boundary stated up front: inferred DM = observer-independent core + observer-relative residual;
we measure the residual, we do NOT claim DM is an artifact.

**Target:** MNRAS. **Carrying figure:** the residual-fraction ladder vs β with the isotropic
control at zero.

---

## 3. PAPER 2 — SIDM anisotropy; needs Phase 1+2 below (cheap)
**Thesis (revised, and now supportable):** *In full 3D, SIDM isotropizes velocity anisotropy 1–2
orders of magnitude faster than the halo collapses, and initial anisotropy produces no
factor-of-two change in gravothermal collapse time; we bound the effect to <X% and show the
collapse is gravothermal (N-independent) rather than relaxation-driven.*

**Do NOT use these two spines — both refuted (§5):** "1D codes can't isotropize" and "fast erasure
undermines the premise."

**Comparison target:** arXiv:2510.23705 (Kamionkowski & Sigurdson), which reports >2× variation
using NSphere-SIDM. Note their own concession of run-to-run variance with "no more than a few
hundred particles in the inner 0.02 kpc" — consistent with our V6 dynamical scatter.

---

## 4. FULL-SCALE AWS PLAN (measured economics)

### 4.1 Cost model — MEASURED, not assumed
- Full step cost (gravity + SIDM) scales as **t ~ N^1.81** [local benchmark, 3 N values]
- Validated against AWS: N=3000 × 4000 steps → model 4.2 min vs observed **3–5 min** compute
  (launch 03:02 → upload 03:08–03:10, minus ~3 min boot+pip measured from the smoke instance)
- **Correction:** high σ/m costs ~1.5–2× more (more scatter pairs): 6000-step σ-ladder cells took
  ~13 min vs model 6.3 min
- **Spot price: $0.0403/hr** for c7i-flex.large (2 vCPU) eu-central-1 → **$0.0202/core-hour**

| N | core-h/run (10⁴ steps) | ×1.5 σ overhead | $/run |
|---|---|---|---|
| 8,000 | 1.0 | 1.6 | **$0.03** |
| 16,000 | 3.6 | 5.5 | $0.11 |
| 30,000 | 11.3 | 17.0 | $0.34 |

### 4.2 THE KEY STRATEGIC FINDING: scale SEEDS, not N
V6 shows scatter does not fall with N over 2000→8000. Therefore **bigger N buys no statistical
power** — it only buys resolution we do not need (V1 already shows the collapse is gravothermal at
N=2000). Power comes from seeds:

seeds needed = (3·s/f)² for a 3σ detection of fractional effect f at scatter s≈12%
→ f=30%: n=2 · **f=10%: n=13** · f=5%: n=52

### 4.3 Phased plan (each phase gated on the previous)
- **Phase 1 — measure the scatter properly (PREREQUISITE).** With n=3 per N, each scatter estimate
  carries ~40% uncertainty; N^−0.21 is NOT distinguishable from N^−0.5. Run **13 seeds × N ∈
  {2000, 8000, 16000}**, arm=none, σ/m=20, 6000 steps = **39 runs ≈ $1.5, ~6 h wall**.
  *Gate:* if scatter falls as ~N^−0.5, revisit large-N; if flat, lock N=8000 forever.
- **Phase 2 — the β campaign at the N Phase 1 selects.** 5 β values × 13 seeds + sham 13 =
  **78 runs ≈ $2.5** at N=8000. Delivers a 10%-level bound on the β-effect.
- **Phase 3 — only if Phase 2 shows a real effect.** σ/m ladder × β at 13 seeds to test whether the
  effect tracks the gravothermal clock. ~130 runs ≈ $4.

**Total for a publishable bound: < $10 and < 24 h wall.** The earlier "130 runs at N=16k, needs a
tree/GPU code" plan was based on an unmeasured cost model and is withdrawn.

### 4.4 Engineering requirements before Phase 2
1. **Checkpoint/resume is mandatory for any run > ~2 h.** A single run is SEQUENTIAL; the fleet
   parallelises across runs, never within one. Spot interruption over a long run loses everything.
   `nbody_fleet.py` already has shard-resume; reuse it rather than re-implement.
2. **Pack 2 cells per instance.** Every cell so far used 1 of 2 vCPUs — free 2× throughput.
3. **Always persist the full ρ_c(t) and β(t) series.** cal1–cal3 did not, and when the collapse
   estimator turned out to be broken those runs could not be re-analysed and had to be re-run.
4. **Smoke-test one instance before every fleet launch.** The first 14-instance launch failed on a
   missing `tqdm`; a 1-instance, ~$0.001 smoke test catches it.
5. Quota: 64 spot + 16 on-demand vCPU = 40 instances. Recently-terminated instances count against
   the quota for several minutes — back off and retry rather than fight `MaxSpotInstanceCountExceeded`.

---

## 5. REFUTED — do not re-propose
1. **"A 1D spherical code cannot isotropize like 3D."** FALSE. arXiv:2506.04334 p.3: NSphere-SIDM
   does a genuine 3D elastic scatter (COM frame, c_α uniform in −1..1, random φ_f), then projects
   the outgoing velocity onto the radial direction and stores c_θ. It isotropizes as 3D does.
2. **"Rapid β erasure undermines their premise."** FALSE. Their setup has t_scat≈110 Myr vs ~16 Gyr
   collapse (ratio 0.007), comparable to ours. β must act as an EARLY IMPRINT on the density/energy
   structure, not as a sustained state.
3. **"L-null is the decisive 3D probe (Panel C)."** FALSE for a spherical halo: plane orientation is
   irrelevant by symmetry, so L-null must be null in 3D too. It is only causal where a preferred
   axis exists (verified: Kerr loss cone, `dmlab kerr_losscone`, asymmetry 0→0.21 as a/M→1).
4. **"N=30,000 is required / the rig needs a tree or GPU code."** FALSE. Based on an unmeasured cost
   model and on assuming Poisson scatter. Measured cost is ~$0.34/run at N=30k, and V6 shows N does
   not reduce scatter anyway.
5. **"t_iso/t_collapse is tightly invariant (0.012–0.014)."** OVERCLAIMED — that came from the
   broken estimator at 3 σ-values. Corrected: 0.005–0.019 (4× spread). Only the order of magnitude
   is supported.

## 6. KNOWN RISKS / OPEN
- Whether NSphere's own β decays as fast as ours is **unverified** — needs their code run or the
  authors' figure. If their β persists longer, that IS a real 1D/3D difference worth chasing.
- All results use one halo (NFW c=10, r_s=0.2) and one σ/m regime; concentration/profile dependence
  untested.
- OM ICs cannot produce β<0; the tangential half of the 1D claim is only reachable via the
  (transient) rotation arm. This is a genuine coverage gap in Paper 2.
