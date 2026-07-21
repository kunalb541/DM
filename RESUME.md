# RESUME — state at laptop restart (2026-07-21 10:53 UTC)

## Everything is safe. The AWS work does NOT depend on this laptop.
Instances run on EC2, write results to S3, and self-terminate. A restart loses nothing
except my local watcher scripts (which only poll and launch — no data lives in them).

## In flight on AWS
- **conv5** (the convergence gate): 4 cells running, checkpointing partial results to
  `s3://nbody-fleet-506250255800-eu-central-1/results/conv5/` every 10 snapshots. ETA was ~10:54 UTC.
  Backstop 330 min, so they cannot outlive it silently.
- **cert1**: 18/24 done, 6 expensive `N8000_s80` cells outstanding →
  `s3://nbody-fleet-506250255800-eu-central-1/results/cert1/`
- **conv4**: DEAD (180-min backstop, 0/18, ~3h lost). Do not wait on it.

## To pick up after restart
```bash
cd ~/Desktop/Research/DM
aws login                                    # session expires; re-auth first
aws s3 ls s3://nbody-fleet-506250255800-eu-central-1/results/conv5/ --region eu-central-1
aws s3 ls s3://nbody-fleet-506250255800-eu-central-1/results/cert1/ --region eu-central-1
aws ec2 describe-instances --region eu-central-1 \
  --filters "Name=tag:Project,Values=sidm-anisotropy" "Name=instance-state-name,Values=running" \
  --query 'length(Reservations[].Instances[])' --output text
```
The filler script (`/tmp/cal/conv5_filler.sh`) dies with the restart and still owed 2 cells;
relaunch those manually if wanted. `/tmp/cal/` is tmp and may be cleared — the launch
artifacts are reproducible from the repo (`fleet/build_package.sh`).

## Blocking state (unchanged)
- Phase 1 BLOCKED until the convergence gate passes.
- All 93 pre-subcycling SIDM results INVALID (see VALIDITY.md); V1–V6 render blank by design.
- Vectorised operator staged, NOT adopted (rejected at 2.02 sigma on beta).
- Open: corrected kappa ≈0.027 > 0.02 criterion → likely lower KAPPA_TARGET to ~0.014.
