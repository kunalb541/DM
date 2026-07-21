# VALIDITY — what may and may not be quoted

**Last updated: 2026-07 (SIDM sub-cycling invalidation).**
Read this before quoting ANY number from `results/` in a paper, plan, or talk.

---

## ⛔ INVALID — all pre-sub-cycling SIDM results (93 files)

Every SIDM result produced before the scattering operator was sub-cycled is **calibration-only
and must not be used for paper claims**. They are tagged in-file with:

```json
"pre_subcycling": true,
"invalid_for_paper": "SIDM scattering under-resolved: kappa ~7x criterion, blocked_frac ~0.52"
```

### Why
A dt-convergence test at fixed physical time (dt = 0.005 / 0.0025 / 0.00125) **failed**:

| symptom | measurement |
|---|---|
| collapse time not converged | t_collapse(K=1.5) moved **−18.4%** from dt=0.005 → 0.0025 (gate: 10–15%) |
| scattering under-resolved | physical scatter rate **26% low** at production dt (28166 → 38218) |
| accuracy criterion violated | kappa = 0.139 / 0.098 / 0.055 = **7.0× / 4.9× / 2.8×** the K&S ≤0.02 |
| saturation bias | `blocked_frac` = **0.52** — the once-per-step cap suppressed half of all candidate scatters |

Root cause: `sidm_scatter` allows each particle at most one scatter per step
(`if scattered[i] or scattered[j]: continue`). At production dt this silently truncates the
scattering rate, non-linearly, and worst exactly where kappa is largest.

### What this invalidates
`V1` gravothermal N-independence · `V2` sham floor · `V3` t_iso/t_collapse ·
`V5` beta-dependence null · `V6` scatter-vs-N. **All require rerunning under sub-cycling.**

---

## ✅ STILL VALID

- **Pure-gravity runs (sigma/m = 0)** — the scattering operator is never invoked, so they are
  untouched. This preserves the `V4` result that equilibrium (Osipkov-Merritt) anisotropy is
  retained **61%** under pure gravity vs **34%** for a rotated IC. (The "5% under SIDM" half of
  V4 is invalidated and must be rerun.)
- **All non-SIDM experiments**: `observer_gamma`, `kerr_losscone`, `stream_shape`, `dm_residual`,
  `causal_locus`, `lyapunov`, `force_law_response`.
- **Method-precision facts** (see README): rotation exactly speed-preserving; KE held to O(1e-4)
  after momentum re-zero; beta-estimator noise ~0.05–0.1 at N ≲ 8000; collapse time must be
  measured by forward threshold crossing.

---

## The fix now in force

`evolve(..., sidm_subcycles='auto')` is the **default**. `sidm_kappa()` probes kappa for a full
step without scattering, then nsub = `ceil(kappa / 0.02)`. A **fixed** nsub is not safe: nsub=8
passed at N=2000 (kappa=0.013) but **failed** at N=3000 (kappa=0.038), because kappa scales with
local density and neighbour count. Runs record `sidm_subcycles_requested`, `sidm_subcycles_used`,
`kappa_full_step`, per-substep `kappa`/`P_max`, and `blocked_frac` (max and mean).

## Rules
1. `analysis.py` **excludes** `pre_subcycling` SIDM by default and prints a loud banner. Do not
   pass `include_pre_subcycling=True` to fill a table.
2. Do not mix pre- and post-sub-cycling SIDM numbers in any figure, table, or sentence.
3. No SIDM claim enters a paper until the `auto` convergence ladder passes AND the headline cells
   have been rerun under sub-cycling.
4. If V1–V6 render blank, that is **correct behaviour**, not a bug.


---

## ⚠️ SECOND DEFECT (2026-07): the recorded `kappa` diagnostic reads ~30% LOW

Found while A/B-testing a vectorised scattering operator. The sequential `sidm_scatter` accumulated
its probability sum **only over pairs it actually evaluated**, skipping pairs blocked by the
once-per-step cap. So the `kappa` written into every result file is an **undercount**.

Measured against an independent probe (`sidm_kappa()`, which sums P over ALL pairs, no scattering):

| source | kappa | ratio to truth |
|---|---|---|
| `sidm_kappa()` (ground truth) | 0.100781 | 1.000 |
| loop diagnostic (what got recorded) | 0.072119 | **0.716** |
| vectorised diagnostic | 0.100781 | 1.000 |

### What this does and does not affect
- **Does NOT affect nsub selection.** `sidm_subcycles='auto'` always called `sidm_kappa()`, which was
  correct. The sub-cycling decisions in every run are sound.
- **DOES affect every recorded `kappa`/`kappa_max` field.** Multiply by **1/0.716 ≈ 1.40** to recover
  the true value.
- **Consequence for the convergence gate:** conv3 reported `kappa_sub = 0.0195` and I called it a
  PASS against the ≤0.02 criterion. Corrected, that is **≈0.027 — it FAILS.** Any convergence claim
  must be re-assessed against corrected kappa, NOT the printed value.
- **The conv4 / cert1 fleets currently running still use the loop operator**, so their `kappa` fields
  will also read ~30% low and need the same correction.

## ⚠️ `blocked_frac` denominator changed
Pre-vectorisation it was (blocked pairs)/(ALL candidate pairs). Post-vectorisation it is
(blocked)/(ACCEPTED pairs) — the meaningful denominator for saturation bias, since only accepted
pairs could have scattered. **The two generations are NOT comparable.** The pre-vectorisation values
quoted elsewhere (0.52 -> 0.066) are on the old denominator.

## Status of the vectorised operator
12.6x faster (17.6 -> 1.4 ms/call) and semantically justified: blocked particles are never
re-evaluated, so pre-computing P from entry velocities reproduces exactly the quantities the
sequential loop would have used. Scatter counts agree (0.27 sigma) and `P_max` is identical.
**NOT YET ADOPTED** - a 5-seed physics check showed beta differing at 2.09 sigma, which is too
marginal to accept or dismiss; a 20-seed test is running. Do not ship it until that passes.
