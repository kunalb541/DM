# CAMPAIGN (pre-registration): Does the anisotropy→SIDM-collapse-time law survive in 3D?

**Status: GATE CLEARED 2026-07. Implementing + validating before fleet launch.**

## 0b. GATE RESULT — ALL THREE CLEAR (verified 2026-07)
- **(a) Paper confirmed.** arXiv:2510.23705, **Kamionkowski & Sigurdson**, "Gravothermal collapse of
  SIDM halos with anisotropic velocity distributions"; collapse times differ from isotropic "by more
  than a factor of two"; uses **Osipkov–Merritt** models — the SAME parameterization our validated
  `nfw_anisotropic_ic.py` produces, so we are directly comparable.
- **(b) NSphere is definitively 1D-spherical** — from the code's own documentation: it reduces "the
  six-dimensional phase space ... to three dimensions (radius r, radial velocity v_r, and angular
  momentum ℓ)", and "the gravitational force on each particle is calculated based solely on the mass
  M(r) enclosed within its current radius (Newton's Shell Theorem)". **Forces are purely radial.**
  The structural argument (cannot host ROI / triaxiality / non-radial torques) is now backed by the
  authors' own method description, not our inference.
- **(c) No 3D study has done this.** Prior gravothermal calibrations were "calibrated to numerical
  simulations initialized with **isotropic** particle velocity distributions"; the leading 3D Arepo
  calibration (arXiv:2504.13004) varies only cross-section, concentration, halo mass — **not**
  initial anisotropy.

### BONUS second contribution (found at the gate)
2510.23705's abstract concedes the effect depends on "global changes in velocity-dispersion profiles
... and **not just on the degree of anisotropy**" — a confound their construction **cannot resolve**,
because their ICs change the whole velocity distribution so β and the σ-profile co-vary. Our
speed-preserving rotation holds each particle's **speed** (hence the global speed/σ magnitude
distribution) fixed and changes only direction — **cleanly separating degree-of-anisotropy from
global-σ-profile change.** This is a second, independent contribution beyond the 1D/3D argument, and
it answers a question the original authors explicitly raise.

### Writeup cautions
- **Naming collision:** arXiv:2504.13004's "β" is the *heat-transfer* parameter, NOT anisotropy β.
  Disambiguate explicitly (use β_aniso vs the conductivity parameter).
- Known OM limitation to acknowledge: OM orbits are "far more radial at large radii than expected in
  cosmological simulations". Our **constant-β rotation arm avoids this** — state that as an advantage.

## 0. VERIFICATION GATE — (retained for the record) do this before spending fleet budget
The entire campaign is framed as a falsification test of **arXiv:2510.23705** ("Gravothermal
collapse of SIDM halos with anisotropic velocity distributions"), which reports t_collapse varying
by >2× (≈3→6 Gyr) as β runs −1/2→+1/2 (radial slower, tangential faster) using **NSphere-SIDM, a 1D
spherically-symmetric shell code**. Our research pass flagged, honestly, that it could **not fully
verify**: the author attribution (came from a search summary, not the arXiv listing), the exact
anisotropy parameterization, and — most important — **whether any 3D N-body study has already done
this**. Absence of evidence ≠ evidence of absence.

**Gate:** dedicated ADS + arXiv full-text search confirming (a) 2510.23705 exists and says what we
claim, (b) its solver really is 1D-spherical, (c) no 3D study already tests SIDM collapse vs initial
anisotropy. If (c) fails, STOP and re-pick. Cheap; do it first.

## 1. The claim we can make that nobody else can
A spherically-symmetric solver **structurally cannot** host the radial-orbit instability (ROI),
triaxial deformation, or non-radial torques. In 3D there are **two competing β-erasing channels** —
ROI (arXiv:1407.6565, 0806.3434) and SIDM collisional isotropization (arXiv:2607.05504; MNRAS 528,
3075) — documented in separate literatures never joined to this question. So: *is the 1D law causal,
attenuated, or reversed in 3D?*

Every published anisotropic result varies β by changing the distribution function, which co-varies
the density realization and virial state — **correlational**. We dial β by **speed-preserving
velocity rotation at fixed positions, and KE / virial ratio fixed to O(1e-4), fixed density profile** — a true
intervention with exact ground truth. That is the rig's monopoly.

## 2. Arms
| arm | what it does | why |
|---|---|---|
| **β-dose (real)** | radialize/tangentialize to span β ∈ [−0.5, +0.5], 7 doses | the causal law |
| **sham** | random-axis rotation, same angle distribution, Δβ≈0 | systematic floor; answers "you just perturbed it out of equilibrium" |
| **L-null** ★ | reorient orbital planes at fixed E, \|L\|, **pericenter** | **THE decisive arm** |
| **N-ladder** | N ∈ {2k,4k,8k,16k}, θ(N)=θ_∞+b·lnN/N on t_collapse | answers "this is two-body relaxation" |
| **σ/m ladder** | sweep cross-section | answers "you inflated σ/m to fit budget" |

**Why L-null is the paper.** It changes 3D orbital-plane geometry (hence ROI susceptibility and
triaxial torque coupling) while leaving *every quantity a 1D spherical code can see* identically
unchanged. **A 1D code predicts EXACTLY ZERO response. A nonzero L-null response is a result no
spherical calculation can produce or explain.** Already validated as exactly-null for spherical
boundaries and causal for orientation-dependent ones (K1, Kerr loss cone).

## 3. Pre-registered kill criteria (locked BEFORE running)
- **KILL the law in 3D:** |Δt_collapse| between β=+0.5 and β=−0.5 < 10% AND consistent with sham at
  2σ after N→∞ extrapolation. (A publishable null — arguably a *stronger* paper than confirmation.)
- **KILL the 3D-specific claim:** L-null response consistent with zero at 2σ → effect is purely
  radial-structure and the 1D treatment is vindicated as sufficient. Report it plainly.
- **Direction is pre-registered:** radial slower, tangential faster (the 1D sign). Recovering the
  **opposite sign in 3D** is the highest-value outcome.
- Report the **differential between matched arms**, never absolute t_collapse.

## 4. Cost
N=16k, ~10⁴ steps. 7 β-doses × 5 seeds × {real, sham} = 70; L-null 4 doses × 5 seeds = 20;
N-ladder 4 × 3 β × 3 seeds = 36. **≈130 production runs.** Budget risk is entirely the 16k tier
(direct O(N²): ~2.6×10⁸ pair-ops/step × 10⁴ steps). N-ladder tiers are cheap. Fleet carries 16k at
reduced seed count if needed. σ/m tuned upward to compress the gravothermal clock into budget —
legitimate and standard, and it directly closes our known SIDM-undersampling gap.

## 5. The carrying figure (three panels, shared β-dose x-axis)
- **A:** t_collapse vs β — the 1D published law (digitized), our 3D real-intervention arm, and the
  sham arm as a shaded systematic band.
- **B:** same, N→∞ extrapolated (relaxation artifacts removed).
- **C:** ★ the L-null arm — a horizontal line at zero labelled "1D prediction," with our measured
  points on it or off it. **Panel C is the paper.**

## 6. Paper
**Thesis:** *The published factor-of-two dependence of SIDM gravothermal collapse time on velocity
anisotropy is derived from a spherically-symmetric solver; under matched-pair causal intervention in
full 3D — where the radial-orbit instability and collisional isotropization both operate — we measure
whether that law is causal, attenuated, or reversed, and we isolate a purely three-dimensional
orbital-plane response that no spherical calculation can produce.*

**Target: MNRAS main journal.** Not PRD (referee pool overlaps the 1D paper's authors → reads as a
methods quibble). Not A&C (risks a third "computational note, not new science" rejection).

**Framing that fixes the two prior rejections:** this is **falsification of a specific 2025 claim**,
NOT "we discovered anisotropy matters." Cite 2510.23705 in the abstract's first sentence. Lead with
the 1D/3D *structural* argument (a spherical code cannot host ROI) — that is a physics argument, not
a resolution argument, and it is unanswerable.

## 7. Reviewer risks → pre-emptions
| risk | pre-emption |
|---|---|
| "Already known (2510.23705)" — **highest risk** | Frame as falsification; structural 1D/3D argument |
| "That's two-body relaxation, not gravothermal" — near-certain at N=16k | θ(N) ladder in **main text**, not appendix; report t_collapse in units of *measured* relaxation time |
| "You inflated σ/m" | σ/m ladder; show β-dependence scales with the gravothermal clock, not σ/m independently; publish the code→physical conversion table |
| "Your rotation makes a non-equilibrium state" | That is the sham arm's entire job; report virial-ratio time series for every run |
| "N=16k too low" | Concede, quantify with the ladder, claim only the **differential** |
| "Idealized, no cosmology" | The claim is about a mechanism **under intervention** — which cosmological sims structurally cannot test. Say once, don't overclaim about real dwarfs. |

## 8. Runs from
`dmlab.py` (unified core: interventions incl. `l_null`, SIDM/feedback/MOND force laws, f_peri,
circularity, β profiles, Jeans, axis ratios) + `nbody_fleet.py` (AWS spot fleet, shard resume,
self-terminating). Add a `sidm_collapse` experiment to the dmlab registry before launch.

---
# CALIBRATION OUTCOME (2026-07, AWS runs cal1-cal5, 73 cells, ~$3) — CAMPAIGN NOT LAUNCHED

## Verdict: the rig cannot currently adjudicate this question. Do NOT run the 130-run campaign.

### Solid results (keep these)
- **Collapse is GRAVOTHERMAL, not two-body relaxation**: t_collapse flat (1150-1350) across N=2000->8000
  while t_relax grows 3.4x; scaling t_collapse ~ N^-0.04. Pre-empts the biggest reviewer objection.
- **Isotropization/collapse RATIO INVARIANCE**: t_iso/t_collapse = 0.013/0.014/0.012 at sigma/m=5/40/80.
  Anisotropy is erased ~75x faster than collapse at ANY cross-section.
- **Equilibrium beta persists under pure gravity (61% at t=1000) but is destroyed by SIDM (5%)** ->
  the erasure channel is collisional isotropization, NOT collisionless phase mixing (sigma/m=0 control).
- Sham floor is SCATTER not bias (0.5 sigma over 7 seeds) - the matched-pair control is sound.

### The null (underpowered)
No beta-dependence of collapse time above sample variance at N=3000, 4 seeds. Critically, the trend's
SIGN FLIPS with the (arbitrary) collapse threshold: slope = +396 (K=1.3), +711 (K=1.5), **-919** (K=2.0),
all ~1 sigma. That sign flip is the signature of noise. A factor-2 effect IS excluded; a <~30% effect is NOT.

### Two hypotheses of mine that were REFUTED (recorded so they aren't re-proposed)
1. **"A 1D code can't isotropize like 3D" - WRONG.** Read arXiv:2506.04334 p.3: NSphere-SIDM performs a
   genuine 3D elastic scatter (COM frame, c_alpha uniform in -1..1, random phi_f) then projects the final
   velocity onto the radial direction and stores c_theta. It isotropizes exactly as 3D does. The
   "dimensional artifact" spine for the paper does not exist.
2. **"Rapid beta erasure undermines their premise" - WRONG.** Their own setup has t_scat~110 Myr vs
   ~16 Gyr collapse (ratio 0.007), comparable to ours. Fast isotropization is common to both, so beta must
   act as an EARLY IMPRINT on the density/energy structure, not as a sustained state.

### The blocker
rho_c is resolved by a handful of particles at N=3000 -> 20-30% scatter and a collapse detector that
FAILED two ways: (a) post-collapse core evacuation puts the global minimum at the END of the run, so the
"first return above initial" search returned NaN for runs that clearly collapsed; (b) marginal runs gave
spurious late crossings (the 3780 "outlier" had max rho_c = 10.36 vs rho_0 = 9.62). Fixed by a forward
threshold-crossing estimator, but the underlying resolution problem remains. Kamionkowski et al. concede
even 100k particles leave "no more than a few hundred particles in the inner 0.02 kpc".
**Required: N >~ 30k and ~10 seeds. At direct O(N^2) that is ~25 h/run on c7i-flex.large -> needs a
tree/FMM or GPU solver before this question is answerable here.**
