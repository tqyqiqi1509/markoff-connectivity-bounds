"""Shared exact checks; old input code and certificates are read-only."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import os
import sys

if sys.flags.optimize:
    raise RuntimeError("Run without -O; certificate checks must remain enabled")
sys.set_int_max_str_digits(0)
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[2]
LEGACY = ROOT / "agents/resumable_exact_determinants"
sys.path.insert(0, str(LEGACY))
import resumable_crt_det as crt
from flint import fmpz_mat, nmod_mat, ctx


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp." + str(os.getpid()))
    temporary.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    temporary.replace(path)


def read_snapshot(path):
    payload = Path(path).read_bytes()
    return json.loads(payload), hashlib.sha256(payload).hexdigest(), payload


def exact_int(value):
    require(type(value) is int or isinstance(value, str), "Non-integer matrix/certificate entry")
    if isinstance(value, str):
        body = value[1:] if value.startswith(("-", "+")) else value
        require(bool(body) and body.isascii() and body.isdecimal(), "Invalid integer string")
    return int(value)


def load_matrix(path, minor_index=0, row_indices=None, column=0, expected_sha=None):
    source, digest, _ = read_snapshot(path)
    if expected_sha:
        require(digest == expected_sha, "Source SHA256 differs from explicitly required bytes")
    raw = source["rows"] if isinstance(source, dict) else source
    require(bool(raw) and bool(raw[0]), "Only nonempty square minors are supported")
    n = len(raw[0])
    require(len(raw) >= n and all(len(row) == n for row in raw), "Invalid source dimensions")
    if row_indices is None:
        require(type(minor_index) is int and minor_index >= 0, "Invalid minor index")
        chosen = list(range(n))
        if minor_index:
            require(len(raw) > n, "No replacement row")
            chosen[(minor_index - 1) % n] = n + (minor_index - 1) % (len(raw) - n)
    else:
        chosen = list(row_indices)
    require(len(chosen) == n and len(set(chosen)) == n
            and all(type(j) is int and 0 <= j < len(raw) for j in chosen), "Invalid ordered row selection")
    require(type(column) is int and 0 <= column < n, "Invalid RHS column")
    rows = [[exact_int(x) for x in raw[j]] for j in chosen]
    binding = {
        "producer_version": 1, "equation": "A Q = D e_i",
        "source_file": str(Path(path)), "source_sha256": digest,
        "shape": [len(raw), n], "row_indices": chosen,
        "column_order": list(range(n)), "column": column,
        "d": source.get("d") if isinstance(source, dict) else None,
    }
    return rows, binding, source


def same_binding(left, right):
    # Paths are diagnostic; ordered bytes, row/column indices and RHS bind the equation.
    return {k: v for k, v in left.items() if k != "source_file"} == {
        k: v for k, v in right.items() if k != "source_file"}


def verify_equation(rows, Q, D, column, pure_python=False):
    n = len(rows)
    require(D > 0 and len(Q) == n, "Invalid positive denominator/vector shape")
    g = D
    for q in Q:
        g = math.gcd(g, q)
    require(g == 1, "Q,D are not primitive")
    if pure_python:
        for j, row in enumerate(rows):
            require(sum(a * q for a, q in zip(row, Q)) == D * int(j == column),
                    "Exact AQ=D e_i failed at row " + str(j))
        engine = "Python arbitrary-precision integer row products"
    else:
        require(fmpz_mat(rows) * fmpz_mat([[q] for q in Q])
                == fmpz_mat([[D * int(j == column)] for j in range(n)]),
                "Exact AQ=D e_i failed")
        engine = "FLINT exact integer matrix multiplication; Python integer gcd"
    return {"all_rows_checked": n, "AQ_equals_D_e_i": True, "positive_primitive_QD": True,
            "arithmetic": engine, "conclusion": "D divides det(A), including possible det(A)=0",
            "nonsingularity_certified": False}


def load_verified_qd(source, checkpoint, pure_python=False):
    record, qd_sha, _ = read_snapshot(checkpoint)
    require(record.get("schema") == "primitive_one_rhs_divisor_v1", "Wrong QD checkpoint schema")
    binding = record["binding"]
    rows, current, src = load_matrix(
        source, row_indices=binding["row_indices"], column=binding["column"],
        expected_sha=binding["source_sha256"])
    require(same_binding(binding, current), "QD source/row/column binding mismatch")
    D = exact_int(record["denominator"])
    Q = [exact_int(x) for x in record["numerators"]]
    verification = verify_equation(rows, Q, D, binding["column"], pure_python)
    return rows, current, src, Q, D, record, qd_sha, verification
