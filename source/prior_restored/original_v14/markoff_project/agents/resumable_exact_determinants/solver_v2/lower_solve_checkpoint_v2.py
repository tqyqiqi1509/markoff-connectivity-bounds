#!/usr/bin/env python3
"""Isolated early-checkpoint stage for the general lower-quotient workflow.

No frozen v1 file is modified. Solve, base determinant, and later lattice
certification are separate phases. A saved candidate is never itself called
a lower-lattice certificate. A legacy solve checkpoint can be imported after
source/row binding and exact BQ=DC checks; its already-computed determinant
is retained with explicit provenance rather than silently recomputed.

Run all expensive phases under the root-owned resource guard. The solver C
call itself is not resumable; only a returned Q,D can be preserved.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from multidim_lower_quotient_v1 import atomic, primes, smooth_strip, fmpz_mat
from flint_den_abi import canonicalize, solve_den, verify_solution
import resumable_crt_det as crt


BINDING_KEYS = ('d', 'source_sha256', 'row_indices', 'restored_content_factors',
                'shape', 'codimension', 'column_order')


def source_lower(path, d):
    payload = path.read_bytes()
    source = json.loads(payload)
    assert source['d'] == d and source['nd'] == 4 * d
    selected = [(i, row, meta) for i, (row, meta) in
                enumerate(zip(source['rows'], source['row_metadata'])) if meta['n'] <= 3 * d]
    factors = [smooth_strip(int(meta['divided_small_content']), primes(6 * d))
               for _, _, meta in selected]
    rows = [[int(x) * f for x in row[:3 * d - 2]]
            for (_, row, _), f in zip(selected, factors)]
    m, n = len(rows), len(rows[0])
    assert n > m and all(len(row) == n for row in rows)
    binding = dict(d=d, source=str(path), source_sha256=hashlib.sha256(payload).hexdigest(),
                   row_indices=[i for i, _, _ in selected],
                   restored_content_factors=factors, shape=[m, n],
                   codimension=n - m, column_order=list(range(n)))
    return [row[:m] for row in rows], [row[m:] for row in rows], binding


def bind(record, binding):
    for key in BINDING_KEYS:
        assert record[key] == binding[key], ('Checkpoint does not bind source', key)


def equation_check(B, C, Q, D, pure_python):
    assert D > 0 and len(Q) == len(B)
    assert all(len(row) == len(C[0]) for row in Q)
    g = D
    for row in Q:
        for x in row:
            g = math.gcd(g, x)
    assert g == 1
    if pure_python:
        return verify_solution(B, C, Q, D)
    assert fmpz_mat(B) * fmpz_mat(Q) == fmpz_mat(C) * D
    return {'BQ_equals_DC': True, 'positive_minimal_common_denominator': True,
            'arithmetic': 'exact FLINT integer matrices and Python integer gcd',
            'nonsingularity': 'separate base determinant required'}


def python_default(B, C):
    Z = fmpz_mat(B).solve(fmpz_mat(C))
    m, r = Z.nrows(), Z.ncols()
    D = 1
    for i in range(m):
        for j in range(r):
            D = math.lcm(D, int(Z[i, j].denominator))
    Q = [[int(Z[i, j].numerator) * (D // int(Z[i, j].denominator))
          for j in range(r)] for i in range(m)]
    return canonicalize(Q, D)


def save_candidate(path, binding, Q, D, timings, provenance):
    record = dict(binding, schema=2, certificate_mode='compressed',
                  status='candidate_QD_saved_equation_not_yet_verified',
                  D=str(D), Q=Q, timings=timings, solution_provenance=provenance)
    # This precedes BQ=DC and every determinant/HNF computation.
    atomic(path, record)
    return record


def verify_complete_crt(matrix_path, state_path, expected_B):
    rows, binding = crt.matrix_input(matrix_path)
    assert rows == expected_B, 'CRT matrix must be exactly this ordered base block'
    reduced, factor, reduction = crt.content_reduce(rows)
    assert factor != 0
    hadamard = crt.hadamard_data(reduced)
    state = json.loads(state_path.read_text())
    assert state['status'] == 'certified_exact'
    residue, modulus = crt.validate_state(state, binding, reduction, hadamard,
                                          recompute=True, A=fmpz_mat(reduced))
    assert int(state['content_factor']) == factor
    assert modulus > 1 << (hadamard['bound_bits'] + 1)
    centered = residue if 2 * residue <= modulus else residue - modulus
    assert abs(centered) <= 1 << hadamard['bound_bits']
    value = factor * centered
    assert int(state['determinant']) == value
    return value, {'method': 'all CRT modular determinants recomputed; proved primes; exact Hadamard bound',
                   'crt_sha256': hashlib.sha256(state_path.read_bytes()).hexdigest()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path)
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--checkpoint', type=Path, required=True)
    ap.add_argument('--phase', choices=['solve', 'verify', 'base-det', 'export-base'], default='solve')
    ap.add_argument('--solver', choices=['python_default', 'dixon_den', 'multi_mod_den'], default='python_default')
    ap.add_argument('--import-legacy', type=Path)
    ap.add_argument('--python-verify', action='store_true')
    ap.add_argument('--base-matrix', type=Path)
    ap.add_argument('--crt-state', type=Path)
    args = ap.parse_args()
    B, C, binding = source_lower(args.source, args.d)
    if args.phase == 'export-base':
        if args.base_matrix is None:
            ap.error('--base-matrix is required for export-base')
        atomic(args.base_matrix, {'rows': B, 'lower_binding': binding,
                                 'purpose': 'exact ordered lower base block for resumable determinant'})
        print(json.dumps({'status': 'base_matrix_exported_no_determinant_claim',
                          'shape': [len(B), len(B)], 'path': str(args.base_matrix)}))
        return

    if args.checkpoint.exists():
        if args.import_legacy:
            ap.error('Do not import over an existing checkpoint; use a new destination')
        rec = json.loads(args.checkpoint.read_text())
        bind(rec, binding)
    elif args.import_legacy:
        old_payload = args.import_legacy.read_bytes()
        old = json.loads(old_payload)
        bind(old, binding)
        Q, D = canonicalize(old['Q'], int(old['D']))
        rec = save_candidate(args.checkpoint, binding, Q, D, {}, {
            'method': 'import prior saved exact solve',
            'legacy_sha256': hashlib.sha256(old_payload).hexdigest(),
            'legacy_file': str(args.import_legacy)})
        if 'base_minor_det' in old:
            value = int(old['base_minor_det'])
            assert value and value % D == 0
            rec.update(base_minor_det=str(value), base_det_provenance={
                'method': 'retained previous exact computation; not independently recomputed on import',
                'legacy_sha256': hashlib.sha256(old_payload).hexdigest()})
    else:
        if args.phase != 'solve':
            ap.error('This phase requires an existing checkpoint or --import-legacy')
        t = time.monotonic()
        if args.solver == 'python_default':
            Q, D = python_default(B, C)
            timings = {'solve_and_export_seconds': time.monotonic() - t}
        else:
            Q, D, timings = solve_den(B, C, args.solver)
        rec = save_candidate(args.checkpoint, binding, Q, D, timings,
                             {'method': args.solver, 'FLINT_version': '3.6.0'})

    Q, D = rec['Q'], int(rec['D'])
    t = time.monotonic()
    rec['solution_verification'] = equation_check(B, C, Q, D, args.python_verify)
    rec.setdefault('timings', {})['latest_equation_check_seconds'] = time.monotonic() - t
    rec['status'] = 'exact_solution_saved_base_det_pending'
    if 'base_minor_det' in rec:
        assert int(rec['base_minor_det']) and int(rec['base_minor_det']) % D == 0
        rec['status'] = 'exact_solution_and_base_det_saved_lattice_not_started'
    atomic(args.checkpoint, rec)

    if args.phase == 'base-det' and 'base_minor_det' not in rec:
        t = time.monotonic()
        if args.crt_state:
            if args.base_matrix is None:
                ap.error('--base-matrix is required with --crt-state')
            value, provenance = verify_complete_crt(args.base_matrix, args.crt_state, B)
        else:
            value = int(fmpz_mat(B).det())
            provenance = {'method': 'exact FLINT fmpz_mat determinant',
                          'resumability': 'Q,D already saved; an interrupted determinant is restarted'}
        assert value != 0 and value % D == 0
        rec['base_minor_det'] = str(value)
        rec['base_det_provenance'] = provenance
        rec['timings']['base_det_phase_seconds'] = time.monotonic() - t
        rec['status'] = 'exact_solution_and_base_det_saved_lattice_not_started'
        atomic(args.checkpoint, rec)
    print(json.dumps({'status': rec['status'], 'shape': binding['shape'],
                      'checkpoint': str(args.checkpoint), 'timings': rec['timings']}), flush=True)


if __name__ == '__main__':
    main()
