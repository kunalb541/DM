# DM — causal N-body experiments on dark-matter halo structure

Idealised, collisionless-primary N-body with **matched-pair causal interventions** and exact ground
truth. The distinguishing instrument is a speed-preserving velocity rotation that holds positions
and (to O(1e-4)) kinetic energy and virial ratio fixed, so orbit structure can be *intervened on*
rather than merely correlated with — something cosmological simulations structurally cannot do.

## Quick start
```bash
python dmlab.py list             # registered experiments
python dmlab.py run <name>       # run one
python dmlab.py verify           # regression: reproduces published numbers
python -m pytest tests/ -q       # method invariants
```

## Layout
| path | contents |
|---|---|
| `dmlab.py` | **the lab** — shared core + experiment registry + CLI. Start here. |
| `nbody_3d.py`, `nbody_dm_ic.py`, `nfw_anisotropic_ic.py`, `sidm.py`, `feedback.py` | physics engine imported by `dmlab` |
| `fleet/` | AWS execution: `nbody_fleet.py`, `cal_cell.py` (one cell per instance), `build_package.sh` |
| `docs/` | **`PAPER_PLAN.md` is the current plan**; plus the idea board and calibration write-ups |
| `results/` | measured fleet output (JSON, with full rho_c(t) and beta(t) series) |
| `tests/` | invariant tests for the causal method |

## `dmlab` core
- **interventions** — `radialize`, `tangentialize`, `sham` (matched null), `l_null` (reorients the
  orbital plane at fixed E, |L| **and pericenter**)
- **measurement** — `gamma` (inner slope), `pericenters`/`f_peri`, `circularity`,
  `beta_sigma_profiles`, `jeans_mass`, `axis_ratios`, `dose_for_beta`
- **force laws** — `newton`, `mond` (algebraic Milgrom), `superfluid` (gated MOND); SIDM scattering
  and potential-fluctuation feedback available in `evolve`

## Running on AWS
```bash
bash fleet/build_package.sh      # -> /tmp/dmlab_pkg.tgz
```
Upload with a cell list, launch spot instances that self-shard on `ami-launch-index`, write JSON to
S3 and self-terminate. **Rules learned the hard way:** smoke-test one instance first (a missing
`tqdm` once killed all 14); always persist the full time series (runs whose series were not saved
could not be re-analysed when an estimator turned out to be wrong); pack 2 cells per instance (each
uses 1 of 2 vCPUs); checkpoint any run over ~2 h against spot interruption.

## Measured cost model
Full step cost scales as `t ~ N^1.81`; spot c7i-flex.large is $0.0403/hr (= $0.0202/core-hour).
A run at N=8000 x 10^4 steps is ~1.6 core-hours ~ $0.03. Run-to-run scatter in collapse time is
**dynamical, not Poisson** (~10-16%, roughly N-independent), so statistical power comes from
**more seeds, not more particles**.

## Method precision (measured — matters for claims)
- The velocity **rotation is exactly speed-preserving** (median |dv|/v = 0 to machine precision).
- The subsequent **momentum re-zero** perturbs individual speeds ~0.3-0.6% and total **KE by ~1e-4**.
  Write "kinetic energy held to O(1e-4)", **not** "fixed by construction". It is common-mode across
  the real and sham arms, so matched-pair *differentials* are unaffected.
- The **beta estimator carries ~0.05-0.1 noise** at N <~ 8000. Quote beta doses with that uncertainty.
- Collapse time must use **forward threshold crossing** (`rho_c >= K*rho_c(0)`), never "first return
  above initial after the global minimum" — post-collapse core evacuation puts the global minimum at
  the end of the run and silently returns NaN for runs that clearly collapsed.
