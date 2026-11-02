#!/usr/bin/env python3
"""Save and verify one exact RHS before any base determinant is available.

Reuses the already-audited solver_v2.python_default; this is not a new solver.
The expensive C solve is not internally resumable. Root must guard large runs.
"""
import argparse
from pathlib import Path
import sys
import time
from common import (LEGACY, atomic, ctx, exact_int, load_matrix, read_snapshot,
                    require, same_binding, sha, stamp, verify_equation)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--checkpoint", type=Path, required=True)
    ap.add_argument("--source-sha256")
    ap.add_argument("--minor-index", type=int, default=0)
    ap.add_argument("--row-indices", type=__import__("json").loads)
    ap.add_argument("--column", type=int, default=0)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--python-verify", action="store_true")
    args = ap.parse_args()
    require(args.threads >= 1, "threads must be positive")
    require(args.checkpoint.resolve() != args.source.resolve(), "Output must not overwrite source")
    ctx.threads = args.threads
    rows, binding, _ = load_matrix(args.source, args.minor_index, args.row_indices,
                                   args.column, args.source_sha256)
    n = len(rows)
    if args.checkpoint.exists():
        record, _, _ = read_snapshot(args.checkpoint)
        require(record.get("schema") == "primitive_one_rhs_divisor_v1", "Wrong checkpoint schema")
        require(same_binding(record["binding"], binding), "Existing QD checkpoint binding mismatch")
        print(__import__("json").dumps({"stage": "reuse_saved_QD", "utc": stamp()}), flush=True)
    else:
        solver_dir = LEGACY / "solver_v2"
        sys.path.insert(0, str(solver_dir))
        from lower_solve_checkpoint_v2 import python_default
        solver_sha = sha(solver_dir / "lower_solve_checkpoint_v2.py")
        start = time.monotonic()
        print(__import__("json").dumps({"stage": "one_RHS_solve_without_base_det",
              "utc": stamp(), "n": n, "threads": args.threads,
              "solver": "existing solver_v2.python_default"}), flush=True)
        Q_matrix, D = python_default(rows, [[int(j == args.column)] for j in range(n)])
        Q = [row[0] for row in Q_matrix]
        record = {
            "schema": "primitive_one_rhs_divisor_v1",
            "status": "candidate_equation_not_yet_verified", "binding": binding,
            "denominator": str(D), "numerators": [str(q) for q in Q],
            "solve_elapsed_seconds": time.monotonic() - start, "created_utc": stamp(),
            "solver": {"entrypoint": "solver_v2.lower_solve_checkpoint_v2.python_default",
                       "source_sha256": solver_sha},
            "base_determinant_required": False, "determinant_computed": False,
        }
        # Save immediately after returned Q,D; verification can be resumed.
        atomic(args.checkpoint, record)
        print(__import__("json").dumps({"stage": "candidate_QD_saved",
              "denominator_bits": D.bit_length(), "utc": stamp()}), flush=True)
    D = exact_int(record["denominator"])
    Q = [exact_int(q) for q in record["numerators"]]
    start = time.monotonic()
    record["verification"] = verify_equation(rows, Q, D, args.column, args.python_verify)
    require(sha(args.source) == binding["source_sha256"], "Source changed during solve")
    record.update(status="verified_exact_divisor_base_determinant_unknown",
                  verified_utc=stamp(), verification_seconds=time.monotonic()-start)
    atomic(args.checkpoint, record)
    print(__import__("json").dumps({
        "status": record["status"], "n": n, "denominator_bits": D.bit_length(),
        "denominator_digits": len(str(D)), "solve_elapsed_seconds": record["solve_elapsed_seconds"],
        "verification_seconds": record["verification_seconds"],
        "checkpoint": str(args.checkpoint), "checkpoint_sha256": sha(args.checkpoint),
        "nonsingularity_certified": False, "determinant_computed": False}), flush=True)


if __name__ == "__main__":
    main()
