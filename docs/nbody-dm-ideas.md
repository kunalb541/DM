# Dark-matter idea board — pre-registered kill criteria

Sharp, falsifiable tests runnable in the existing rig. **Kill criteria are locked here
BEFORE running** (pre-registration): if a test hits its kill line, the idea is dropped,
recorded honestly. Ordered cheapest-to-kill first. Unique instruments available:

- Matched-pair, speed-preserving velocity rotations (positions fixed; KE/Q to O(1e-4)).
- The `f_peri` handle (fraction with pericenter < r_c) — validated causal knob on central mass.
- Three coring engines in one rig: orbit-reshaping, SIDM (`sidm.py`), feedback (`feedback.py`).
- Collisionality law θ(N)=θ_∞+b·lnN/N; merger virial-vs-escape lead/lag.

| # | Claim | Kill criterion | Cost | Status |
|---|-------|----------------|------|--------|
| 1 | Apparent inner slope γ depends on tracer orbit type (observer-dependent core-cusp) | \|Δγ(radial−tangential)\| < 0.1 → drop | zero sims | **DONE** |
| 2 | f_peri is the universal handle mediating ALL coring (orbit, SIDM, feedback) | SIDM cores at FIXED f_peri → dead | 1 run + collapse | **DONE — FALSIFIED as universal (→ #5)** |
| 3 | Merger virial-vs-escape lead/lag changes under SIDM (timing lag as a DM probe) | lead/lag unchanged across σ/m → drop | hours | todo |
| 4 | Radial-orbit-instability is a real collisionless coring channel for cores/anisotropic systems | undershoot relaxes back OR vanishes as N→∞ → drop | matched sweep | todo |
| 5 | SIDM vs feedback cores separable by CAUSAL response sign (not just inner slope) | same response sign/magnitude → drop | 2 cores + kick | todo |
| 6 | f_peri × σ/m coupling: radial halos core faster under SIDM at fixed σ/m | SIDM core-time independent of initial f_peri → drop | matched pairs + SIDM | todo |
| 7 | Predictive≠causal inversion at gravothermal collapse (density predicts, T-gradient causes) | same observable tops both rankings → drop | ~10⁴-step runs | todo |
| 8 | Collisionality thermometer: place real cored dwarfs on θ_∞+b·lnN/N to read granularity | real halos don't lie on the curve → drop | (longshot) | todo |

## Batch B — radical, rig-testable (pre-registered 2026-07)

| # | Claim | Kill criterion | Cost | Status |
|---|-------|----------------|------|--------|
| B1 | The cusp is causally NON-LOCAL — outer plunging orbits (f_peri pop) set the inner slope, not inner residents | outer-shell (r>2r_s) radialize moves γ <20% as much per-particle as inner shell → LOCAL | 5 short runs → 30 (rigorous) | **DONE — non-local NOT supported (local-to-r_s)** |
| B2 | Halo hysteresis — identical present-day (ρ,σ,β) haloes have different futures (non-Markovian; DM halo = a history not a state) | two histories → same endpoint evolve identically within sham noise → Markovian | matched pair | todo |
| B3 | Coring has a measurable irreversibility timescale (Loschmidt echo: v→−v after evolving) | return fidelity stays high for all T → reversible, no entropy barrier | cheap | **DONE — chaotic, mean-field (λ_∞=0.34)** |
| B4 | MOND vs DM distinguishable by CAUSAL response to orbit reshaping (breaks the rotation-curve fit degeneracy) | identical Δγ/Δf_peri response signatures → degeneracy is causal too | modify force kernel | **DONE — NULL (its "hint" refuted by K4)** |
| B5 | Inferred SHAPE (triaxiality c/a) is observer-dependent by tracer orbit type (geometric analogue of orig-#1 γ) | \|Δ(c/a)(radial−tangential)\| < 0.05 → shape observer-invariant | 1 ROI run | **DONE — PASS (Δc/a=0.44)** |
| B6 | "Best predictor" is TASK-relative — no globally-privileged summary statistic | one variable tops every target's ranking (γ, outer slope, c/a, T_relax) → privileged descriptor exists | multi-target rank | todo |

## Batch K — literature-driven, rig-testable (pre-registered 2026-07)

| # | Claim | Kill criterion | Cost | Status |
|---|-------|----------------|------|--------|
| K1 | L-null (pericenter-preserving plane rotation) stays NULL for spherical sink but becomes CAUSAL for a Kerr (spin-oriented) loss cone | χ=0 asymmetry not ~0, OR χ→1 asymmetry <0.01 → no break | semi-analytic | **DONE — BREAKS (A: 0→0.21)** |
| K2 | SIDM gravothermal collapse is a predictive≠causal INVERSION: density predicts, T-gradient/heat-flux causes | same observable tops predictive AND interventional rankings → no inversion | ~10⁴-step SIDM | todo |
| K3 | FDM granule heating (correlated potential noise) has a DIFFERENT chaos spectrum than CDM lnN/N granularity (coherence-length/time diffusion law) | FDM λ-scaling matches CDM lnN/N form → no distinct fingerprint | noise model + Loschmidt | todo |
| K4 | Superfluid/MOND-like force changes WHICH intervention direction is causal (asymmetry flip), not the slope | response-asymmetry sign same as Newtonian across gating → no switch | force kernel + gate | **DONE — NULL (0.1σ; refutes B4's hint)** |
| K5 | Stream-inferred halo shape depends on the stream's orbit class, not the halo alone | radial vs tangential stream tracers infer same c/a & slope → no bias | stream tracers, fixed halo | **DONE — PASS (b/a 0.71→0.90)** |
| K6 | f_peri is the orbit-space shadow of phase-space folding; fold count/parity beats density for future central response | density/veldisp top predictive+causal over fold count/parity → no | fold-topology predictors | todo |
| K7 | A ring caustic produces observer-dependent false cusp/core or stream-gap WITHOUT a bound subhalo | ring caustic leaves inferred slope/shape unchanged → no lensing artifact | inject ring potential | todo |

## Results log (honest, as they land)

- **#1 — DONE (observer-dependent γ).** SURVIVES, but confound-corrected: a naive angular-
  momentum split gave Δγ=0.85 (inflated — L correlates with current radius); the clean
  **circularity** split (η=L/L_circ(E), radius-independent) gives **Δγ = +0.24** (radial
  tracers γ≈1.96 cuspier, tangential γ≈1.72 flatter). Real, non-negligible observer-dependence
  (= the mass-anisotropy degeneracy that limits dwarf-DM measurements), but a **modest ~0.24
  systematic, NOT a cusp↔core flip** — the strong "core-cusp is an observer artifact" is an
  overclaim. Script: `dmlab run observer_gamma`.

- **#2 — DONE (survives kill, but only co-movement shown).** SIDM (σ/m=40, N=5000, 600
  steps) does **NOT** core at fixed f_peri: as γ drops 1.62(CDM)→1.25(SIDM), f_peri(<0.05)
  drops 0.085→0.060 (~30%) and f_peri(<0.1) drops 0.235→0.206. So scattering isotropizes
  velocities → raises pericenters → drains the low-pericenter population *while* it cores.
  The pre-registered kill line (cores at FIXED f_peri) is not hit, and the **direction matches
  the orbit handle** (lower f_peri ⇒ de-concentrated). **BUT** this is *co-movement, not
  mediation*: it does not prove f_peri *causes* the coring vs. being a correlated by-product.
  The strong "universal handle" claim needs the collapse test — do orbit + SIDM + feedback all
  fall on ONE Δγ-vs-Δf_peri curve? That's the honest next step. Script: `dmlab run fingerprints`.

- **#2 collapse test — DONE, universal handle FALSIFIED (this becomes result #5).** Put all three
  engines on the Δγ-vs-Δf_peri plane (orbit interventions trace the ground-truth handle curve).
  **They do NOT collapse onto one curve** (the script's auto-"COLLAPSE" verdict was a fitting
  artifact — a free-intercept line + inflated σ absorbed vertical excursions; caught on double-check):
  - **orbit** radialize slope≈1.0, tangentialize slope≈4.5 — f_peri moves γ (its own mechanism), but
    already nonlinear/asymmetric.
  - **SIDM** cores at *nearly fixed* f_peri: Δf_peri≈−0.01…−0.06 but Δγ≈−0.37…−0.48 (needs Δf_peri≈−0.4
    to explain via the handle; overshoots the orbit curve 1.4–6.8×). Large **direct** channel.
  - **feedback** cores with the *wrong sign*: Δf_peri≈+0.01 (up) while Δγ≈−0.08…−0.36 (down). Cannot
    act through f_peri at all.
  **Conclusion:** f_peri is the handle for the *orbit* channel, NOT a universal mediator. The three
  engines have **distinct phase-space fingerprints** (= idea #5: mechanisms are separable by causal
  Δγ/Δf_peri signature, not just inner slope). SIDM/feedback carry coring physics invisible to a
  pericenter-only description — ties to the reinterpretation program (velocity-space info a coarse
  observer misses). Also corrects the earlier #2 "survives" read: SIDM's f_peri shift is real but far
  too small to *mediate* its coring. Scripts: `dmlab run fingerprints`, fig `outputs/figures/collapse_test.png`.

- **B5 — DONE, PASS (`dmlab run shape_observer`).** Observer-dependent SHAPE confirmed and large. Drove a
  radial-orbit-unstable NFW (r_a=0.10) triaxial (global c/a 0.98→0.55 over 1500 steps), then split by
  orbit type (labels assigned at t=0 while spherical): radial tracers c/a=0.48 (elongated bar),
  tangential tracers c/a=0.92 (nearly round), **Δ(c/a)=−0.44** (~9× the kill). Generalizes the
  observer-dependence from density (orig-#1 γ) to geometry. Honest: expected given ROI physics (radial
  orbits ARE the bar) — a strong confirmation-with-mechanism, not a shock; tangential subset small (5%).
- **B1 — DONE, non-local NOT supported (`dmlab run causal_locus`).** Rigorous rerun (per-shell matched sham +
  3 IC seeds) KILLED the exciting first-pass result. Δγ/1000 by shell: [0.05,0.10) −0.06±0.29 (the
  pass-1 "flattening sign-flip" was NOISE — 301 particles), [0.10,0.20) +0.27±0.07, [0.20,0.40)
  +0.28±0.08, [0.40,0.80) **+0.06±0.06 (consistent with 0)**. Causal control of γ is concentrated at
  r≈0.1–0.4 (around r_s); the fully-outside shell has no significant effect → **cusp control is
  local-to-r_s, NOT non-local.** Pass-1's "non-local + sign flip" was an unmatched-control + single-
  realization artifact — a clean double-check catch (would have reported new physics and been wrong).

- **B3 — DONE (`dmlab run lyapunov`, `dmlab run lyapunov`).** Perturbed cusp is CHAOTIC (two-trajectory
  divergence, clean exponential, R²=0.99). λ(N)=0.549/0.462/0.398 at N=2000/4000/8000 → θ(N)-style fit
  λ=0.340+55.7·(lnN/N), R²=0.994 → **λ_∞=0.34>0: mean-field chaos SURVIVES the collisionless limit**
  (t_lyap≈11.8 t_cross), with a granular *enhancement* at finite N (not a two-body artifact — pure
  two-body would drop λ 0.30× over the range; observed 0.72×). Irreversibility is real but SLOW
  (t_lyap~12 t_cross) and only LOGARITHMIC in observer resolution: t_mix(δ)~ln(r_s/δ)·t_lyap, so
  reversibility is exponentially cheap to buy back by refining δ. **Four double-check catches in B3
  alone:** (1) v1 machine-precision reversal gave F=1 trivially (roundoff seed too small); (2) v2 had
  a broken δ=0.1>eps column; (3) N=16000 direct run too slow → lean λ(N) scan; (4) script's crude
  "GRANULAR→vanishes" verdict WRONG — proper N→∞ extrapolation shows λ_∞>0 (mean-field). ODD reading:
  the arrow is real and observer-resolution-relative. Ties to the θ(N) collisionality program.
- **B4 — DONE (`dmlab run force_law_response`).** MOND (Milgrom algebraic, a0=1) vs Newtonian causal response to
  radialize/tangentialize (each vs its own sham): Newton +0.529/−0.395, MOND +0.443/−0.471. Response
  magnitudes differ only ~0.08 (<0.10 kill) → **NOT a clean discriminator at this strength.** Two real
  hints: MOND baseline cusp steeper (γ 1.69 vs 1.51, expected — MOND boosts gravity), and the response
  ASYMMETRY FLIPS (Newton radialize-dominant, MOND tangentialize-dominant). No error bars → hint, not
  result. Caveats: algebraic (QUMOND-lite) MOND + matched-pair response on a shared IC (not MOND-native
  equilibrium). Worth a seeded follow-up before any claim.

- **K1 — DONE, BREAKS as predicted (`dmlab run kerr_losscone`).** The L-null intervention (reorient orbital
  plane at FIXED E, |L|, pericenter) is EXACTLY null for a spherical/Schwarzschild sink (capture
  asymmetry retro−pro = 0.0000 at χ=0) but becomes causal for a Kerr sink, with a clean monotonic
  dose-response in spin: A(χ) = 0.044, 0.093, 0.159, **0.208** at χ=0.3, 0.6, 0.9, 0.998 (retrograde
  captured more — larger Kerr marginally-bound radius). Same pericenter-preserving intervention that is
  null for spherical is a sensitive spin-orientation probe — the causal toolkit detects frame-dragging
  invisible to the pericenter/density observables. Semi-analytic (orbit-level loss cone; no dynamical
  refilling); r_mb(μ) interpolated linearly between exact equatorial prograde/retrograde Kerr values;
  built from textbook Kerr geodesics, NOT the unverified arXiv:2606.18050. Extends the predictive≠causal
  thread into a GR boundary.

- **K5 — DONE, PASS (`dmlab run stream_shape`).** In a fixed triaxial log potential (1:0.9:0.7), each orbit
  family's inferred in-plane flattening b/a reports the plane IT samples, not the halo's 3D shape:
  x-y tube→0.88(≈qy 0.90), y-z tube→0.76, box→0.71(≈qz 0.70), outer→0.90(rounder). Spread 0.18 → a
  single stream infers its orbit-class-biased shape. **Double-check catch:** my first headline metric
  (footprint c/a) was confounded by orbital planarity (all orbits thin → c/a 0.15-0.26, not the halo
  shape); the clean estimator is b/a (the in-plane ratio). Confirms the claim for the right reason.
  Honesty: footprint proxy for full stream-track fitting; qualitative family-dependence is what's shown.

- **K4 — DONE, NULL (`dmlab run force_law_response`); RETROACTIVELY KILLS B4's hint.** Three force laws
  (Newton / MOND-global / superfluid-gated MOND inside r_sf=0.5), same radialize/tangentialize/sham
  battery, **3 IC seeds for error bars** (which B4 lacked). Asymmetry A=|radialize|−|tangential|:
  Newton +0.171±0.158, MOND-global −0.014±0.187, superfluid +0.139±0.207. **Every A is consistent with
  zero, and Newton vs MOND-global differ by only 0.76σ** → B4's apparent "response-asymmetry flip"
  (Newton radialize-dominant → MOND tangentialize-dominant) was NOISE. |A(Newton)−A(superfluid)|=0.032
  (0.1σ) → no phase switch. **Conclusion: orbit-reshaping causal response does NOT discriminate
  Newton/MOND/superfluid at this strength.** Scope caveat: a0=1.0 vs halo accels ~1–3 = only mildly
  MOND; a deeper-MOND setup (higher a0 / lower-density halo) is the honest follow-up before calling the
  discriminator dead in general. Lesson: flagging B4 as "hint, needs seeds" was right — the seeds refuted it.

## Notes / references
Kerr loss-cone orientation (arXiv:2606.18050, unverified); SIDM gravothermal collapse (arXiv:2506.06269);
FDM granule heating (arXiv:2412.09908); superfluid DM (Berezhiani-Khoury arXiv:1507.03013); X-Stream
(arXiv:2508.02666); cosmic origami (Neyrinck MNRAS 2012); caustic rings (arXiv:1508.04494).
core-cusp a matter of perspective (arXiv:1707.06303); gravothermal collapse robust to feedback;
resonant-SIDM core collapse; FDM dynamical heating. The on-brand thread: dark matter as the next
system for the predictive-vs-causal / observer-dependent program.
