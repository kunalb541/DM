#!/usr/bin/env python3
"""
nbody_fleet.py — S3-backed, spot-safe, reproducible fleet harness
================================================================================
Runs the collisionality-crossover N-ladder study (PREREG_collisionality_
crossover.md) across many AWS spot / on-demand instances that may be killed at
any time. Wraps the REAL kernels — nbody_stress.run_stress (predictive track) and
nbody_aws_battery._cell_worker (causal track) — no physics reinvented.

Design goals (this is exactly what the harness guarantees):
  • SAVE TO BUCKETS  — results stream to an S3 bucket (or a local dir for
    testing); one immutable object per shard. `s3://bucket/prefix` -> S3;
    any other path -> local filesystem. Same code path either way.
  • REPRODUCIBLE     — every seed is a deterministic function of the cell, so a
    re-run reproduces identical rows. A run-level manifest pins git commit, a
    content hash of the kernel source, package versions, and the full grid spec;
    resuming validates the spec hash so a bucket can never mix incompatible runs.
  • EASY TO PICK & PARALLELISE — launch M instances, give each
    `--shard-index k --shard-count M`; each grabs a disjoint 1/M slice of the
    pending shards with NO coordination service. Add instances any time.
  • SPOT-KILL SAFE   — work unit is a small shard (default 25 seeds). A shard is
    done only when its `.done` marker is written (last, after the rows). A kill
    mid-shard leaves rows but no marker; on resume the shard reloads its rows,
    recomputes only the missing seeds, and finishes. Bounded loss = one shard,
    and rows are flushed to the bucket every `--flush-every` seeds so even that
    is small. Default mode is --resume.

Subcommands:
  plan       print the grid, per-cell reps, per-N time config, projected cost
  run        execute this instance's slice of pending shards
  list       show pending / done shard counts (per cell)
  aggregate  concatenate all shard rows for a run into one local jsonl (+analyse)

Examples:
  # one box, everything, local dir (smoke)
  python nbody_fleet.py run --bucket ./outputs/fleet_local --run-id dev --smoke
  # 8 spot instances, S3, DM + stylised, full predictive ladder:
  python nbody_fleet.py run --bucket s3://my-bucket/nbody --run-id cc1 \
      --track predictive --inits nfw3d hernquist3d plummer3d \
      --shard-index $K --shard-count 8 --workers 30
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
import math
import os
import socket
import subprocess
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

# ── Grid definition (keep in sync with PREREG_collisionality_crossover.md) ──────
LADDER: List[int] = [512, 1024, 2048, 4096, 8192, 16384]

# Predictive replicate taper and causal matched-pair taper.
REP_TAPER: Dict[int, int]    = {512: 500, 1024: 500, 2048: 500,
                                4096: 200, 8192: 100, 16384: 50}
CAUSAL_PAIRS: Dict[int, int] = {512: 100, 1024: 100, 2048: 100,
                                4096: 50,  8192: 50,  16384: 50}

# Per-N integration refinement to hold the energy-conservation gate (<1e-2).
# Calibrated by nbody_calibrate_nladder.py: dt=0.005/steps=600 passes for N<=4096
# (drift<=4.2e-3) but FAILS at N=8192 (2.2e-2). STEP_MULT multiplies the step
# count and divides dt (and shifts the early/mid snapshot indices) so the PHYSICAL
# horizon t=steps*dt is identical across the ladder. 16384 value is provisional
# pending the steps=1200 dt-probe; set to 2 if steps=1200 clears the gate there.
STEP_MULT: Dict[int, int] = {512: 1, 1024: 1, 2048: 1, 4096: 1, 8192: 2, 16384: 4}
_BASE_STEPS, _BASE_DT, _BASE_HE, _BASE_HM = 600, 0.005, 100, 300

DEFAULT_INITS  = ["nfw3d", "hernquist3d", "plummer3d"]   # DM first
DEFAULT_MODELS = ["direct_isolated", "pm_periodic"]
DEFAULT_EPS    = [0.05]
BOX_SIZE, PLUMMER_A, NFW_C, K_FINE = 2.0, 0.20, 10.0, 16
CAUSAL_SAFE_NMAX = 4096   # causal worker uses fixed dt=0.005 (see note in run())

# Kernel source files whose content is hashed into the reproducibility manifest.
_KERNEL_SRC = ["nbody_3d.py", "nbody_stress.py", "nbody_dm_ic.py",
               "nbody_aws_battery.py", "nbody_fleet.py"]


# ── time config: hold physical horizon fixed while refining dt at high N ────────
def time_config(n: int) -> Dict:
    m = STEP_MULT.get(n, 1)
    return {"steps": _BASE_STEPS * m, "dt": _BASE_DT / m,
            "h_early": _BASE_HE * m, "h_mid": _BASE_HM * m, "step_mult": m}


# ── storage backends ────────────────────────────────────────────────────────────
class LocalBackend:
    def __init__(self, root: str):
        self.root = os.path.abspath(root)
        os.makedirs(self.root, exist_ok=True)

    def _p(self, key: str) -> str:
        return os.path.join(self.root, key)

    def exists(self, key: str) -> bool:
        return os.path.exists(self._p(key))

    def get_text(self, key: str) -> Optional[str]:
        p = self._p(key)
        if not os.path.exists(p):
            return None
        with open(p, "r") as f:
            return f.read()

    def put_text(self, key: str, text: str) -> None:
        p = self._p(key)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        tmp = p + ".tmp"
        with open(tmp, "w") as f:
            f.write(text)
        os.replace(tmp, p)   # atomic

    def list_keys(self, prefix: str) -> List[str]:
        base = self._p(prefix)
        out: List[str] = []
        if not os.path.exists(base):
            return out
        for dirpath, _dirs, files in os.walk(base):
            for fn in files:
                full = os.path.join(dirpath, fn)
                out.append(os.path.relpath(full, self.root))
        return out

    def describe(self) -> str:
        return f"local:{self.root}"


class S3Backend:
    def __init__(self, bucket: str, prefix: str):
        try:
            import boto3  # lazy: only needed for s3:// URIs
        except ImportError as e:
            raise SystemExit("boto3 required for s3:// buckets — `pip install boto3`") from e
        self.bucket = bucket
        self.prefix = prefix.strip("/")
        self._s3 = boto3.client("s3")

    def _k(self, key: str) -> str:
        return f"{self.prefix}/{key}" if self.prefix else key

    def exists(self, key: str) -> bool:
        from botocore.exceptions import ClientError
        try:
            self._s3.head_object(Bucket=self.bucket, Key=self._k(key))
            return True
        except ClientError:
            return False

    def get_text(self, key: str) -> Optional[str]:
        from botocore.exceptions import ClientError
        try:
            r = self._s3.get_object(Bucket=self.bucket, Key=self._k(key))
            return r["Body"].read().decode("utf-8")
        except ClientError:
            return None

    def put_text(self, key: str, text: str) -> None:
        self._s3.put_object(Bucket=self.bucket, Key=self._k(key),
                            Body=text.encode("utf-8"))

    def list_keys(self, prefix: str) -> List[str]:
        paginator = self._s3.get_paginator("list_objects_v2")
        full_prefix = self._k(prefix)
        strip = len(self.prefix) + 1 if self.prefix else 0
        out: List[str] = []
        for page in paginator.paginate(Bucket=self.bucket, Prefix=full_prefix):
            for obj in page.get("Contents", []):
                out.append(obj["Key"][strip:])
        return out

    def describe(self) -> str:
        return f"s3://{self.bucket}/{self.prefix}"


def make_backend(uri: str):
    if uri.startswith("s3://"):
        rest = uri[len("s3://"):]
        bucket, _, prefix = rest.partition("/")
        return S3Backend(bucket, prefix)
    return LocalBackend(uri)


# ── grid / cells / shards ────────────────────────────────────────────────────────
def _safe(v) -> str:
    return str(v).replace(".", "p").replace("/", "_")


def cell_id(cell: Dict) -> str:
    return (f"{cell['init']}__{cell['model']}__N{cell['n']}"
            f"__eps{_safe(cell['eps'])}")


def build_cells(track: str, inits: List[str], models: List[str],
                eps_list: List[float], ladder: List[int]) -> List[Dict]:
    cells: List[Dict] = []
    if track == "causal":
        models = ["direct_isolated"]   # causal battery is direct-isolated only
    for init in inits:
        for model in models:
            for n in ladder:
                for eps in eps_list:
                    reps = (REP_TAPER if track == "predictive" else CAUSAL_PAIRS).get(n)
                    if not reps:
                        continue
                    cells.append({"track": track, "init": init, "model": model,
                                  "n": n, "eps": float(eps), "reps": reps})
    return cells


def cell_seed_base(cell: Dict) -> int:
    """Deterministic, cell-specific seed base -> reproducible & cross-cell-disjoint."""
    key = f"{cell['track']}|{cell['init']}|{cell['model']}|{cell['n']}|{cell['eps']}"
    h = hashlib.blake2b(key.encode(), digest_size=6).digest()
    # keep well away from other batteries' seed ranges (they use small ints / 2000+)
    return 1_000_000 + int.from_bytes(h, "big") % 500_000_000


def cell_shards(cell: Dict, shard_size: int) -> List[Tuple[int, int]]:
    """List of (seed_start_offset, count) shards covering the cell's reps."""
    reps = cell["reps"]
    return [(s, min(shard_size, reps - s)) for s in range(0, reps, shard_size)]


def shard_prefix(run_id: str, cell: Dict) -> str:
    return f"{run_id}/{cell['track']}/{cell_id(cell)}"


def shard_keys(run_id: str, cell: Dict, start: int, count: int) -> Tuple[str, str]:
    base = f"{shard_prefix(run_id, cell)}/shard_{start:06d}_{count:03d}"
    return base + ".jsonl", base + ".done"


# ── picklable per-seed workers (top-level for ProcessPoolExecutor) ──────────────
def _run_predictive(task: Tuple) -> Dict:
    from nbody_stress import StressConfig, run_stress
    init, model, n, eps, seed, steps, dt, he, hm = task
    cfg = StressConfig(model=model, init=init, seed=int(seed), n=int(n),
                       steps=int(steps), dt=float(dt), eps=float(eps),
                       box_size=BOX_SIZE, k_fine=K_FINE, h_early=int(he),
                       h_mid=int(hm), plummer_a=PLUMMER_A, nfw_c=NFW_C)
    row = run_stress(cfg, use_numba=True)
    row["seed"] = int(seed)
    row["dt"] = float(dt)          # record the per-N refined timestep for provenance
    return row


def _run_causal(task: Tuple) -> Dict:
    from nbody_aws_battery import _cell_worker
    init, _model, n, eps, seed, *_ = task
    row = _cell_worker((init, int(n), float(eps), int(seed)))
    return row


# ── shard execution with resume ─────────────────────────────────────────────────
def _load_existing_rows(backend, jsonl_key: str) -> Tuple[List[Dict], set]:
    txt = backend.get_text(jsonl_key)
    if not txt:
        return [], set()
    rows, seeds = [], set()
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue   # tolerate a torn final line from a mid-write kill
        rows.append(r)
        if "seed" in r:
            seeds.add(int(r["seed"]))
    return rows, seeds


def run_one_shard(backend, run_id: str, cell: Dict, start: int, count: int,
                  workers: int, flush_every: int, meta: Dict) -> str:
    jsonl_key, done_key = shard_keys(run_id, cell, start, count)
    if backend.exists(done_key):
        return "skip"

    tc = time_config(cell["n"])
    base = cell_seed_base(cell)
    target_seeds = [base + start + i for i in range(count)]

    rows, have = _load_existing_rows(backend, jsonl_key)
    todo = [s for s in target_seeds if s not in have]
    if not todo:
        _write_done(backend, done_key, cell, start, count, len(rows), meta)
        return "resume-complete"

    worker = _run_predictive if cell["track"] == "predictive" else _run_causal
    tasks = [(cell["init"], cell["model"], cell["n"], cell["eps"], s,
              tc["steps"], tc["dt"], tc["h_early"], tc["h_mid"]) for s in todo]

    from nbody_3d import _worker_init
    done_since_flush = 0
    with ProcessPoolExecutor(max_workers=workers, initializer=_worker_init,
                             initargs=(True,)) as ex:
        futs = {ex.submit(worker, t): t for t in tasks}
        for fut in as_completed(futs):
            try:
                rows.append(fut.result())
            except Exception:
                rows.append({"seed": futs[fut][4], "status": "error",
                             "message": traceback.format_exc()})
            done_since_flush += 1
            if done_since_flush >= flush_every:
                backend.put_text(jsonl_key, "\n".join(json.dumps(r) for r in rows) + "\n")
                done_since_flush = 0

    backend.put_text(jsonl_key, "\n".join(json.dumps(r) for r in rows) + "\n")
    _write_done(backend, done_key, cell, start, count, len(rows), meta)
    return "done"


def _write_done(backend, done_key: str, cell: Dict, start: int, count: int,
                n_rows: int, meta: Dict) -> None:
    marker = {"cell_id": cell_id(cell), "track": cell["track"], "n": cell["n"],
              "start": start, "count": count, "n_rows": n_rows,
              "time_config": time_config(cell["n"]),
              "seed_base": cell_seed_base(cell),
              "host": meta.get("host"), "git": meta.get("git"),
              "finished_utc": _utcnow()}
    backend.put_text(done_key, json.dumps(marker, indent=2))


# ── reproducibility manifest ─────────────────────────────────────────────────────
def _utcnow() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__))).decode().strip()
    except Exception:
        return "unknown"


def _code_hash() -> str:
    h = hashlib.blake2b(digest_size=16)
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in sorted(_KERNEL_SRC):
        p = os.path.join(here, fn)
        if os.path.exists(p):
            with open(p, "rb") as f:
                h.update(f.read())
    return h.hexdigest()


def _pkg_versions() -> Dict:
    out = {"python": sys.version.split()[0]}
    for pkg in ("numpy", "scipy", "numba"):
        try:
            out[pkg] = __import__(pkg).__version__
        except Exception:
            out[pkg] = "absent"
    return out


def grid_spec(track: str, inits, models, eps_list, ladder) -> Dict:
    return {"track": track, "inits": sorted(inits), "models": sorted(models),
            "eps": sorted(eps_list), "ladder": sorted(ladder),
            "rep_taper": REP_TAPER, "causal_pairs": CAUSAL_PAIRS,
            "step_mult": STEP_MULT, "box_size": BOX_SIZE,
            "plummer_a": PLUMMER_A, "nfw_c": NFW_C, "k_fine": K_FINE}


def _spec_hash(spec: Dict) -> str:
    return hashlib.blake2b(json.dumps(spec, sort_keys=True).encode(),
                           digest_size=12).hexdigest()


def ensure_manifest(backend, run_id: str, spec: Dict, force_new: bool) -> Dict:
    key = f"{run_id}/manifest.json"
    existing = backend.get_text(key)
    spec_h = _spec_hash(spec)
    if existing and not force_new:
        man = json.loads(existing)
        if man.get("spec_hash") != spec_h:
            raise SystemExit(
                f"Grid spec mismatch for run '{run_id}':\n"
                f"  bucket spec_hash = {man.get('spec_hash')}\n"
                f"  this-run  spec_hash = {spec_h}\n"
                "Refusing to mix incompatible results in one bucket. Use a new "
                "--run-id, or --force-new-manifest to overwrite (loses linkage).")
        return man
    man = {"run_id": run_id, "created_utc": _utcnow(), "git": _git_commit(),
           "code_hash": _code_hash(), "packages": _pkg_versions(),
           "spec": spec, "spec_hash": spec_h}
    backend.put_text(key, json.dumps(man, indent=2))
    return man


# ── pending-shard enumeration + modular slicing ─────────────────────────────────
def all_shards(cells: List[Dict], shard_size: int) -> List[Tuple[Dict, int, int]]:
    out: List[Tuple[Dict, int, int]] = []
    for cell in cells:
        for (start, count) in cell_shards(cell, shard_size):
            out.append((cell, start, count))
    # deterministic order so modular slicing is stable across instances
    out.sort(key=lambda t: (cell_id(t[0]), t[1]))
    return out


def pending_shards(backend, run_id: str, shards) -> List[Tuple[Dict, int, int]]:
    pend = []
    for (cell, start, count) in shards:
        _jsonl, done_key = shard_keys(run_id, cell, start, count)
        if not backend.exists(done_key):
            pend.append((cell, start, count))
    return pend


# ── subcommands ──────────────────────────────────────────────────────────────────
def cmd_plan(args) -> None:
    cells = build_cells(args.track, args.inits, args.models, args.eps, args.ladder)
    calib = None
    if os.path.exists(args.calibration):
        calib = json.load(open(args.calibration))
    print(f"track={args.track}  cells={len(cells)}  shard_size={args.shard_size}")
    print(f"{'cell_id':<52}{'reps':>6}{'steps':>7}{'dt':>9}{'shards':>8}")
    total_sim = 0
    per_sim_s = (calib or {}).get("per_sim_s", {})
    est_core_h = 0.0
    for cell in cells:
        tc = time_config(cell["n"])
        ns = len(cell_shards(cell, args.shard_size))
        total_sim += cell["reps"]
        print(f"{cell_id(cell):<52}{cell['reps']:>6}{tc['steps']:>7}"
              f"{tc['dt']:>9.4f}{ns:>8}")
        t1 = per_sim_s.get(cell["model"], {}).get(str(cell["n"]))
        if t1:
            mult = 2.0 if cell["track"] == "causal" else 1.0
            est_core_h += cell["reps"] * t1 * tc["step_mult"] * mult / 3600.0
    print(f"\ntotal sims/pairs = {total_sim}")
    if est_core_h:
        print(f"projected ≈ {est_core_h:.1f} core-h "
              f"(from {args.calibration}; step_mult-scaled)")
        for v in (24, 48, 96):
            print(f"   on {v:>3} vCPU ≈ {est_core_h / v:.2f} h wall")
    else:
        print("(run nbody_calibrate_nladder.py to populate a cost projection)")


def cmd_list(args) -> None:
    backend = make_backend(args.bucket)
    cells = build_cells(args.track, args.inits, args.models, args.eps, args.ladder)
    shards = all_shards(cells, args.shard_size)
    pend = pending_shards(backend, args.run_id, shards)
    print(f"backend={backend.describe()}  run={args.run_id}")
    print(f"shards: total={len(shards)}  done={len(shards) - len(pend)}  "
          f"pending={len(pend)}")
    by_cell: Dict[str, List[int]] = {}
    for (cell, _s, _c) in pend:
        by_cell.setdefault(cell_id(cell), []).append(1)
    for cid, lst in sorted(by_cell.items()):
        print(f"  pending {len(lst):>4}  {cid}")


def cmd_run(args) -> None:
    backend = make_backend(args.bucket)
    if args.smoke:
        args.inits = args.inits[:1]
        args.models = ["direct_isolated"]
        args.ladder = [min(args.ladder)]
        args.shard_size = min(args.shard_size, 3)
        # cap reps so smoke is a handful of sims, not the full taper
        global REP_TAPER, CAUSAL_PAIRS
        REP_TAPER = {n: min(v, 6) for n, v in REP_TAPER.items()}
        CAUSAL_PAIRS = {n: min(v, 6) for n, v in CAUSAL_PAIRS.items()}

    spec = grid_spec(args.track, args.inits, args.models, args.eps, args.ladder)
    man = ensure_manifest(backend, args.run_id, spec, args.force_new_manifest)
    meta = {"host": socket.gethostname(), "git": man["git"]}

    cells = build_cells(args.track, args.inits, args.models, args.eps, args.ladder)
    if args.max_reps:
        # Cap reps for this invocation (quick tests / budget runs). Seeds are a
        # deterministic prefix, so a later full run extends the same range — the
        # capped rows stay valid and are reused on resume.
        for c in cells:
            c["reps"] = min(c["reps"], args.max_reps)
    if args.track == "causal":
        unsafe = [c for c in cells if c["n"] > CAUSAL_SAFE_NMAX]
        if unsafe and not args.force_causal_highn:
            print(f"NOTE: dropping {len(unsafe)} causal cells with N>{CAUSAL_SAFE_NMAX} "
                  "— the causal worker uses fixed dt=0.005, which fails the "
                  "conservation gate at high N. Extend with a dt-calibrated causal "
                  "worker, then pass --force-causal-highn. (Predictive track has no "
                  "such limit — its dt is refined per-N via STEP_MULT.)")
            cells = [c for c in cells if c["n"] <= CAUSAL_SAFE_NMAX]

    shards = all_shards(cells, args.shard_size)
    pend = pending_shards(backend, args.run_id, shards)
    # disjoint 1/M slice for this instance
    mine = [s for i, s in enumerate(pend) if i % args.shard_count == args.shard_index]

    print(f"backend={backend.describe()}  run={args.run_id}  host={meta['host']}")
    print(f"instance {args.shard_index+1}/{args.shard_count}: "
          f"{len(mine)} of {len(pend)} pending shards ({len(shards)} total)")
    if args.smoke:
        print("SMOKE MODE: 1 init, direct_isolated, smallest N, tiny shards")

    done = skipped = 0
    for (cell, start, count) in mine:
        status = run_one_shard(backend, args.run_id, cell, start, count,
                               args.workers, args.flush_every, meta)
        if status == "done":
            done += 1
            print(f"  [done]  {cell_id(cell)} shard {start:+d}×{count}")
        else:
            skipped += 1
    print(f"finished: {done} shards computed, {skipped} already complete. "
          f"Re-run with --resume (default) after any spot kill to continue.")


def cmd_aggregate(args) -> None:
    backend = make_backend(args.bucket)
    cells = build_cells(args.track, args.inits, args.models, args.eps, args.ladder)
    shards = all_shards(cells, args.shard_size)
    out_rows: List[Dict] = []
    n_missing = 0
    for (cell, start, count) in shards:
        jsonl_key, done_key = shard_keys(args.run_id, cell, start, count)
        if not backend.exists(done_key):
            n_missing += 1
            continue
        rows, _ = _load_existing_rows(backend, jsonl_key)
        out_rows.extend(rows)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        for r in out_rows:
            f.write(json.dumps(r) + "\n")
    print(f"aggregated {len(out_rows)} rows from {len(shards)-n_missing}/{len(shards)} "
          f"complete shards -> {args.out}"
          + (f"  ({n_missing} shards still pending)" if n_missing else ""))
    if args.analyse and args.track == "predictive":
        from nbody_stress import analyse
        res = analyse(out_rows, n_boot=1000)
        ana_path = args.out.replace(".jsonl", "_analysis.json")
        json.dump({k: v for k, v in res.items()}, open(ana_path, "w"),
                  indent=2, default=str)
        print(f"wrote analysis -> {ana_path}")


# ── CLI ──────────────────────────────────────────────────────────────────────────
def _add_common(p) -> None:
    p.add_argument("--track", choices=["predictive", "causal"], default="predictive")
    p.add_argument("--inits", nargs="+", default=DEFAULT_INITS)
    p.add_argument("--models", nargs="+", default=DEFAULT_MODELS)
    p.add_argument("--eps", type=float, nargs="+", default=DEFAULT_EPS)
    p.add_argument("--ladder", type=int, nargs="+", default=LADDER)
    p.add_argument("--shard-size", type=int, default=25)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_plan = sub.add_parser("plan", help="show grid + projected cost")
    _add_common(p_plan)
    p_plan.add_argument("--calibration", default="outputs/nladder_calibration_steps600.json")
    p_plan.set_defaults(func=cmd_plan)

    p_list = sub.add_parser("list", help="show pending/done shards")
    _add_common(p_list)
    p_list.add_argument("--bucket", required=True)
    p_list.add_argument("--run-id", required=True)
    p_list.set_defaults(func=cmd_list)

    p_run = sub.add_parser("run", help="run this instance's slice of shards")
    _add_common(p_run)
    p_run.add_argument("--bucket", required=True,
                       help="s3://bucket/prefix  OR  a local directory path")
    p_run.add_argument("--run-id", required=True)
    p_run.add_argument("--shard-index", type=int, default=0)
    p_run.add_argument("--shard-count", type=int, default=1)
    p_run.add_argument("--workers", type=int,
                       default=max(1, (os.cpu_count() or 2) - 1))
    p_run.add_argument("--flush-every", type=int, default=10)
    p_run.add_argument("--max-reps", type=int, default=0,
                       help="cap reps per cell for this invocation (0 = full taper)")
    p_run.add_argument("--smoke", action="store_true")
    p_run.add_argument("--force-new-manifest", action="store_true")
    p_run.add_argument("--force-causal-highn", action="store_true")
    p_run.set_defaults(func=cmd_run)

    p_agg = sub.add_parser("aggregate", help="concatenate shard rows -> one jsonl")
    _add_common(p_agg)
    p_agg.add_argument("--bucket", required=True)
    p_agg.add_argument("--run-id", required=True)
    p_agg.add_argument("--out", default="outputs/fleet_aggregate.jsonl")
    p_agg.add_argument("--analyse", action="store_true")
    p_agg.set_defaults(func=cmd_aggregate)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
