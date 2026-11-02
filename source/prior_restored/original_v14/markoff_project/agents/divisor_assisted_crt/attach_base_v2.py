#!/usr/bin/env python3
"""Attach a separately certified nonzero base determinant to saved Q,D.

No solve/determinant is rerun. Output exactly matches the current one-cofactor
v2 binding. A base file's mathematical certification is a separate premise.
"""
import argparse
import hashlib
import json
from pathlib import Path
from common import atomic, load_verified_qd, read_snapshot, require, sha, stamp


def attach(source, qd_file, base_file, output, expected_base_sha=None, pure_python=False):
    require(not Path(output).exists(), "Write a new v2 destination; preserve previous input bytes")
    require(Path(output).resolve() not in {Path(source).resolve(), Path(qd_file).resolve(), Path(base_file).resolve()},
            "Output must not overwrite an input")
    rows, b, src, Q, D, original, qd_sha, verification = load_verified_qd(
        source, qd_file, pure_python)
    n = len(rows)
    require(b["row_indices"] == list(range(n)) and b["column_order"] == list(range(n)),
            "Existing one_cofactor v2 supports only the first n source rows and original column order")
    require(isinstance(src, dict) and "d" in src and len(src["rows"]) > n,
            "Existing one_cofactor v2 requires the original overdetermined source and d")
    base, base_sha, _ = read_snapshot(base_file)
    if expected_base_sha:
        require(base_sha == expected_base_sha, "Frozen base certificate SHA mismatch")
    require(base["matrix_sha256"] == b["source_sha256"], "Base source SHA mismatch")
    entry = base["minor_certificates"][0]
    require(entry["row_indices"] == list(range(n)), "Wrong ordered base minor")
    determinant = int(entry["det"])
    require(determinant != 0 and determinant % D == 0, "Need separately certified nonzero base determinant divisible by D")
    # Match code/one_cofactor_column_certificate.py literally, including path text.
    binding = {
        "method_version": 2, "equation": "A Q = D e_i",
        "source_sha256": b["source_sha256"],
        "base_minors_sha256_at_start": base_sha,
        "base_determinant_sha256": hashlib.sha256(str(determinant).encode()).hexdigest(),
        "source_file": str(source), "n": n, "column": b["column"], "d": src["d"],
    }
    record = {
        "status": "exact_cofactor_equation_verified", "binding": binding,
        "denominator": str(D), "numerators": [str(q) for q in Q],
        "solve_elapsed_seconds": original["solve_elapsed_seconds"],
        "created_utc": stamp(), "verified_utc": stamp(),
        "import_provenance": {
            "original_QD_file": str(qd_file), "original_QD_sha256": qd_sha,
            "no_new_solve_or_determinant": True,
            "base_determinant": "Imported from separately certified frozen base file",
            "equation_rechecked": verification,
            "source_path_text": "Use exactly this source argument in the existing v2 producer"},
    }
    require(sha(source) == b["source_sha256"] and sha(base_file) == base_sha
            and sha(qd_file) == qd_sha, "An input changed during v2 attachment")
    atomic(output, record)
    return {"status": record["status"], "output": str(output),
            "output_sha256": sha(output), "base_sha256": base_sha}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--qd", type=Path, required=True)
    ap.add_argument("--base-minors", type=Path, required=True)
    ap.add_argument("--base-sha256")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--python-verify", action="store_true")
    a = ap.parse_args()
    print(json.dumps(attach(a.source, a.qd, a.base_minors, a.output,
                            a.base_sha256, a.python_verify)), flush=True)


if __name__ == "__main__":
    main()
