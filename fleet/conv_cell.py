#!/usr/bin/env python3
"""SIDM dt-convergence cell. Same PHYSICAL time at each dt; records safety diagnostics.
Usage: conv_cell.py <tag> <N> <dt> <sigma> <seed> <phys_steps_at_dt0>"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, dmlab as D

tag, N, dt, sig, seed, base_steps = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), \
    float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
_n = sys.argv[7] if len(sys.argv) > 7 else 'auto'
NSUB = _n if _n == 'auto' else int(_n)
steps = int(base_steps * (0.005/dt))          # fixed physical time
snap  = max(int(20 * (0.005/dt)), 1)          # fixed physical cadence
cfg = D.cfg_for(N, dt=dt); mass = 1.0/N
pos, vel = D.sample(N, seed, kind='isotropic')
vel = D.radialize(pos, vel, 0.245)            # beta ~ +0.5
b0 = D.mean_beta(pos, vel)
rng = np.random.default_rng(seed*13+5); diag = {}
p, v = pos.copy(), vel.copy()
ts, rcs, bets = [], [], []
import subprocess
BUCKET = os.environ.get('DMLAB_BUCKET', 'nbody-fleet-506250255800-eu-central-1')
RUNID  = os.environ.get('DMLAB_RUN', 'conv5')
def _partial(ts, rc, bet, done):
    """Upload progress so far. Without this, a backstop shutdown loses the ENTIRE run --
    which is exactly what killed conv4 after ~3h of compute."""
    rec = dict(tag=tag, N=N, dt=dt, sigma=sig, seed=seed, steps_done=done, steps_target=steps,
               partial=True, sidm_subcycles_used=diag.get('sidm_subcycles_used'),
               kappa_full_step=diag.get('kappa_full_step'), P_max=diag.get('P_max'),
               kappa_max=diag.get('kappa'), blocked_frac=diag.get('blocked_frac'),
               n_scatter_total=diag.get('n_scatter_total'), beta_i=b0,
               t_series=ts, rho_series=rc, beta_series=bet)
    open(f'/tmp/{tag}.partial.json','w').write(json.dumps(rec))
    subprocess.run(['aws','s3','cp',f'/tmp/{tag}.partial.json',
                    f's3://{BUCKET}/results/{RUNID}/{tag}.partial.json','--region','eu-central-1',
                    '--only-show-errors'], check=False)

for k in range(0, steps, snap):
    p, v = D.evolve(p, v, cfg, mass, snap, 'newton', dt=dt, sigma_over_m=sig, rng=rng, sidm_diag=diag,
                    sidm_subcycles=NSUB)
    r = np.linalg.norm(p - D.CEN, axis=1)
    ts.append((k+snap)*dt)                     # PHYSICAL time, comparable across dt
    rcs.append(float((r < 0.12).sum())/N/(4/3*np.pi*0.12**3))
    bets.append(D.mean_beta(p, v))
    if (k//snap) % 10 == 9:                     # checkpoint to S3 every 10 snapshots
        _partial([(kk+snap)*dt for kk in range(0, k+1, snap)], rcs, bets, k+snap)
ts = np.array(ts); rc = np.array(rcs); bet = np.array(bets)
def tcol(K):
    j = np.where(rc >= K*rc[0])[0]
    return float(ts[j[0]]) if len(j) else float('nan')
t_iso = float('nan')
if abs(b0) > 0.05:
    j = np.where(np.abs(bet) <= 0.5*abs(b0))[0]
    if len(j): t_iso = float(ts[j[0]])
rec = dict(tag=tag, N=N, dt=dt, sigma=sig, seed=seed, steps=steps, beta_i=b0,
           sidm_subcycles_requested=NSUB, sidm_subcycles_used=diag.get('sidm_subcycles_used'),
           kappa_full_step=diag.get('kappa_full_step'),
           blocked_frac_mean=diag.get('blocked_frac_mean'), dt_eff=dt/max(diag.get('sidm_subcycles_used') or 1, 1),
           t_collapse_K13=tcol(1.3), t_collapse_K15=tcol(1.5), t_collapse_K20=tcol(2.0),
           t_iso=t_iso, rho_max=float(rc.max()), rho_final=float(rc[-1]),
           scatters_per_step=diag.get('n_scatter_total', 0)/max(steps, 1),
           P_max=diag.get('P_max'), kappa_max=diag.get('kappa'),
           blocked_frac=diag.get('blocked_frac'), n_pairs=diag.get('n_pairs'),
           t_series=ts.tolist(), rho_series=rc.tolist(), beta_series=bet.tolist())
os.makedirs('/tmp/cal', exist_ok=True)
open(f"/tmp/cal/{tag}.json", "w").write(json.dumps(rec))
print("CELL_DONE", tag, "kappa=%.4f P_max=%.3f blocked=%.3f" %
      (rec['kappa_max'] or 0, rec['P_max'] or 0, rec['blocked_frac'] or 0), flush=True)
