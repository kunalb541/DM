> **STATUS: FUTURE WORK — NOT THE CURRENT PLAN.**
> The canonical plan is [`PLAN_papers_and_scaleup.md`](PLAN_papers_and_scaleup.md).
> This paper's spine (R3 / figure F3, the Δf_peri–Δγ fingerprints) is **under-seeded**: one
> realization per point, no error bars. Do not submit or frame around it until the seeded
> battery (`dmlab run fingerprints`) is complete.

# THE PAPER — one paper, planned from scratch

## Title
**Causal fingerprints of dark-matter core formation: cores are not made by rearranging orbits,
and the mechanisms that make them are separable in phase space**

## One-sentence thesis
Using matched-pair, speed-preserving velocity interventions with exact ground truth, we show that a
cored profile is **not dynamically accessible** from an NFW cusp by redistributing orbits alone, and
that the three mechanisms which *can* produce cores — orbit reshaping, self-interaction, and
repeated potential fluctuation — act through **causally distinct channels** that occupy separable
loci in the (Δf_peri, Δγ) plane, giving each a signature that survives the degeneracy of the final
density profile.

## Why this is a dark-matter paper (not a methods note)
The core–cusp problem and the **diversity problem** (halos of the same mass with very different
inner densities) both hinge on one unanswered question: *which mechanism made a given core?* Final
profiles are degenerate — SIDM, feedback, and orbit structure can all end at similar γ. We supply
what the final profile cannot: **how each mechanism moves a halo through phase space.** Two halos
with identical γ are shown to be distinguishable by the causal route that produced them.

---

## The argument in five results

| # | result | status | data |
|---|--------|--------|------|
| R1 | **Accessibility null:** orbit rearrangement alone cannot core an NFW cusp — ±1–2% central-mass modulation, bidirectional, dose-ordered, sham-null, surviving N→∞ extrapolation (θ_∞ = +12.7±0.7 / −9.7±0.6 ×10⁻³) | **DONE** | fleet `cc1` |
| R2 | **The handle:** within the intervention class f_peri is necessary *and* sufficient for the central-mass response (sufficiency battery, multivariate transfer, pericenter-preserving null ≤0.5%) | **DONE** | fleets `ft2`,`cc1` |
| R3 | ★ **f_peri is NOT universal — the spine.** Three engines occupy distinct loci: orbit moves *along* the f_peri axis (slope≈1); **SIDM cores at ~fixed f_peri** (Δf_peri≈−0.03 vs Δγ≈−0.4, off the orbit curve by 1.4–6.9×); **feedback cores with wrong-sign f_peri** (raises it while γ falls) | **NEEDS SEEDS** | `dmlab run fingerprints` (1 realization/point) |
| R4 | **Why: the mechanism behind SIDM's fingerprint.** SIDM erases anisotropy *collisionally* — β retained 61% under pure gravity vs **5%** under scattering — at 50–200× the collapse rate. It therefore flattens γ by isotropising velocities **without draining low-pericenter orbits**, which is exactly the fixed-f_peri signature of R3 | **DONE** | `cal4`,`cal6` |
| R5 | **Controls.** SIDM coring is gravothermal, not two-body relaxation: d(ln t_col)/d(ln N) = +0.01/−0.08/−0.06 at K=1.3/1.5/2.0 over N=2000→8000 while t_relax grows 3.4×. Sham floor is scatter not bias (0.2–1.6σ) | **DONE** | `cal6` |

**R4 is what lifts R3 from an empirical curiosity to a mechanism.** SIDM's fixed-f_peri fingerprint
is *predicted* by its isotropisation channel, and we measure both.

---

## Figures

| fig | content | data status |
|-----|---------|-------------|
| **F1** | Accessibility null: central-mass response vs dose for radialize / tangentialize / sham, with the N→∞ extrapolation inset | ready |
| **F2** | f_peri as causal handle: response vs Δf_peri, dose-ordered, sham at zero | ready |
| **F3** ★ | **THE FIGURE.** (Δf_peri, Δγ) plane: orbit arm traces the handle curve through the origin; SIDM points sit vertically off it; feedback points sit on the wrong side. Error bars from seeds | **needs the battery** |
| **F4** | β(t) with and without scattering (61% vs 5%) + t_iso/t_collapse across σ/m — the mechanism behind F3 | ready |
| **F5** | Controls: t_collapse vs N (flat) against t_relax (rising 3.4×) | ready |

Only **F3** requires new compute.

---

## The one battery still needed (≈$5, <1 day)

Power comes from **seeds, not particles** — measured: run-to-run scatter is dynamical
(sd/mean ≈ 10–16%, ~N^−0.21, not Poisson N^−0.5), so large N buys nothing. Cost is $0.03/run at
N=8000 (measured t~N^1.81, spot $0.0403/hr).

```
orbit arm     radialize θ ∈ {0.2,0.4,0.8}, tangentialize θ ∈ {0.2,0.4,0.8}   × 8 seeds =  48
SIDM arm      σ/m ∈ {10,20,40,80}                                            × 8 seeds =  32
feedback arm  amp ∈ {0.1,0.2,0.4}                                            × 8 seeds =  24
sham          matched to each arm's rotation magnitude                       × 8 seeds =  24
N-ladder      N ∈ {2k,4k,8k,16k} on the three fingerprints                   × 3 seeds =  36
                                                                             TOTAL ≈ 164 runs
```
Each cell reports (Δγ, Δf_peri) against its matched sham, plus the full ρ_c(t)/β(t) series.
Engineering rules already learned: persist all series; smoke-test one instance first; pack 2 cells
per instance; checkpoint any run >2 h.

---

## Target and framing

**MNRAS main journal.** Not A&C (twice burned by "computational note"), not PRD.

**Framing that pre-empts the two prior rejections:**
- Lead with **accessibility** ("can a core be *reached*?") and **mechanism separability**, NOT with
  "anisotropy affects the inner slope" — that is Hansen / Del Popolo and is exactly what drew the
  "well-known facts" rejection.
- The novelty is **causal**: interventions with exact ground truth, which no cosmological simulation
  can perform. We report matched-pair *differentials*, never absolute profiles.
- Disclose AI assistance explicitly (the declaration text already drafted).

**Anticipated referee objections → pre-emption**
| objection | answer |
|---|---|
| "Anisotropy–slope link is known" | We claim accessibility and separability, not the link. R1 is a null; R3 is new. |
| "Your coring is two-body relaxation" | R5: N-independent while t_relax grows 3.4×, at three thresholds. |
| "The intervention makes a non-equilibrium state" | That is the sham arm's entire job: floor is 0.2–1.6σ (R5). |
| "σ/m is inflated" | Conceded explicitly. R3 claims **signature separability**, not physical core-formation rates. |
| "Idealised isolated halo" | The claim is about a mechanism *under intervention*, which cosmological sims structurally cannot test. Stated once, no overclaiming about real dwarfs. |

## Honest limitations (in the paper, not buried)
- Feedback is a toy potential-fluctuation model; two variants bracket the truth.
- SIDM/feedback run at inflated σ/m and short timescales — separability claim only.
- One halo (NFW c=10); concentration dependence untested.
- OM ICs cannot produce β<0.

## Explicitly NOT in this paper
Observer-dependence suite (Δγ=0.24, Δc/a=0.44, stream bias, Jeans residual 37–184%), the Kerr
L-null result, the SIDM-anisotropy null, and the reinterpretation program. All stay in the repo as
future work. One paper.
