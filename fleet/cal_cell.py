#!/usr/bin/env python3
"""One calibration cell -> its own output file (for parallel xargs execution)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import dmlab as D
tag, N, steps, sig, arm, dose, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), \
    float(sys.argv[4]), sys.argv[5], float(sys.argv[6]), int(sys.argv[7])
r = D.exp_sidm_collapse(N=N, steps=steps, sigma_over_m=sig, arm=arm, dose=dose,
                        seeds=(seed,), snap=20)
run = r['runs'][0]
rec = {k: run[k] for k in ('seed','beta_i','beta_f','t_iso','t_collapse','rho_min','t_min','ca_final','t_col_k')}
rec['beta_series']=run['beta_series']; rec['t_series']=run['t_series']; rec['rho_series']=run['rho_series']
rec.update(tag=tag, N=N, steps=steps, sigma=sig, arm=arm, dose=dose)
open(f"/tmp/cal/{tag}.json","w").write(json.dumps(rec))
print("CELL_DONE", tag, flush=True)
