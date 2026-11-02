#!/usr/bin/env python3
"""Pinned FLINT 3.6.0 public integer common-denominator solve interface.

This owns separate C allocations. It never reads a Python extension object's
layout. The public *_den solvers still use fmpq_mat internally in FLINT 3.6.0;
this is an explicit-algorithm interface, not a rational-allocation bypass.
"""
from __future__ import annotations
import ctypes as C
import math
from pathlib import Path
import sys
import time

if sys.flags.optimize:
    raise RuntimeError('Run without -O: exact certificate checks must remain enabled')

HERE = Path(__file__).resolve().parent
DEPS = HERE.parents[1] / 'computational_mechanisms' / 'deps'
sys.path.insert(0, str(DEPS))
import flint  # Load its pinned GMP/MPFR dependency chain before ctypes.


class FmpzMat(C.Structure):
    # Official v3.6.0 src/fmpz_types.h, not a Python extension layout.
    _fields_ = [('entries', C.POINTER(C.c_long)), ('r', C.c_long),
                ('c', C.c_long), ('stride', C.c_long)]


def _library():
    assert C.sizeof(C.c_long) == C.sizeof(C.c_void_p) == 8
    assert flint.__version__ == '0.9.0', 'ABI audited only for python-flint 0.9.0'
    candidates = list((DEPS / 'python_flint.libs').glob('libflint-*.so.*'))
    if len(candidates) != 1:
        raise RuntimeError('Expected exactly the audited vendored FLINT library')
    lib = C.CDLL(str(candidates[0]))
    # flint_version is char[], NOT char*. Never use c_char_p.in_dll here.
    first = C.c_char.in_dll(lib, 'flint_version')
    version = C.string_at(C.addressof(first)).decode('ascii')
    assert version == '3.6.0', ('Unaudited FLINT ABI', version)
    pm, pz = C.POINTER(FmpzMat), C.POINTER(C.c_long)
    signatures = {
        'fmpz_mat_init': ([pm, C.c_long, C.c_long], None),
        'fmpz_mat_clear': ([pm], None),
        'fmpz_mat_entry': ([pm, C.c_long, C.c_long], pz),
        'fmpz_init': ([pz], None),
        'fmpz_clear': ([pz], None),
        'fmpz_set_str': ([pz, C.c_char_p, C.c_int], C.c_int),
        'fmpz_sizeinbase': ([pz, C.c_int], C.c_size_t),
        'fmpz_get_str': ([C.c_void_p, C.c_int, pz], C.c_void_p),
        'fmpz_mat_solve_dixon_den': ([pm, pz, pm, pm], C.c_int),
        'fmpz_mat_solve_multi_mod_den': ([pm, pz, pm, pm], C.c_int),
    }
    for name, (args, result) in signatures.items():
        fn = getattr(lib, name)
        fn.argtypes, fn.restype = args, result
    return lib


def canonicalize(Q, D):
    """Return the positive least common denominator for this rational array."""
    if D == 0:
        raise ValueError('A zero denominator is not a solution')
    if D < 0:
        Q, D = [[-x for x in row] for row in Q], -D
    g = D
    for row in Q:
        for x in row:
            g = math.gcd(g, x)
    if g != 1:
        Q, D = [[x // g for x in row] for row in Q], D // g
    return Q, D


def solve_den(B, C_rhs, algorithm='dixon_den'):
    """Solve B Q = D C exactly via a pinned public C API; caller verifies.

    A nonzero return from FLINT is only a finder result. Certificate callers
    must check BQ=DC and certify nonsingularity separately. Cleanup happens
    even after a Python exception. A killed C call cannot yield a checkpoint.
    """
    if algorithm not in ('dixon_den', 'multi_mod_den'):
        raise ValueError(algorithm)
    m = len(B)
    if not m or any(len(row) != m for row in B):
        raise ValueError('B must be nonempty and square')
    r = len(C_rhs[0])
    if len(C_rhs) != m or not r or any(len(row) != r for row in C_rhs):
        raise ValueError('Incompatible nonempty right hand side')
    lib = _library()
    owned = []
    denominator = C.c_long()
    lib.fmpz_init(C.byref(denominator))
    started = time.monotonic()

    def new_matrix(rows, cols):
        obj = FmpzMat()
        lib.fmpz_mat_init(C.byref(obj), rows, cols)
        owned.append(obj)
        return obj

    def set_matrix(obj, rows):
        for i, row in enumerate(rows):
            for j, value in enumerate(row):
                ptr = lib.fmpz_mat_entry(C.byref(obj), i, j)
                if lib.fmpz_set_str(ptr, str(value).encode('ascii'), 10):
                    raise ValueError('FLINT rejected a decimal integer')

    def integer(ptr):
        # Caller-owned buffer avoids crossing allocator ownership boundaries.
        length = lib.fmpz_sizeinbase(ptr, 10) + 4
        buf = C.create_string_buffer(length)
        if not lib.fmpz_get_str(buf, 10, ptr):
            raise RuntimeError('fmpz_get_str failed')
        return int(buf.value)

    try:
        A, R, X = new_matrix(m, m), new_matrix(m, r), new_matrix(m, r)
        set_matrix(A, B)
        set_matrix(R, C_rhs)
        setup_seconds = time.monotonic() - started
        t = time.monotonic()
        success = getattr(lib, 'fmpz_mat_solve_' + algorithm)(
            C.byref(X), C.byref(denominator), C.byref(A), C.byref(R))
        solve_seconds = time.monotonic() - t
        if not success:
            raise ArithmeticError('FLINT solver did not return a solution')
        t = time.monotonic()
        D = integer(C.byref(denominator))
        Q = [[integer(lib.fmpz_mat_entry(C.byref(X), i, j))
              for j in range(r)] for i in range(m)]
        export_seconds = time.monotonic() - t
        t = time.monotonic()
        Q, D = canonicalize(Q, D)
        normalize_seconds = time.monotonic() - t
        return Q, D, dict(setup_seconds=setup_seconds, solve_seconds=solve_seconds,
                         export_seconds=export_seconds, normalize_seconds=normalize_seconds,
                         total_solver_seconds=time.monotonic() - started)
    finally:
        for obj in reversed(owned):
            lib.fmpz_mat_clear(C.byref(obj))
        lib.fmpz_clear(C.byref(denominator))


def verify_solution(B, C_rhs, Q, D):
    """Independent Python integer verification, including least positive D."""
    m, r = len(B), len(C_rhs[0])
    assert D > 0 and len(Q) == m and all(len(row) == r for row in Q)
    columns = list(zip(*Q))
    for i, row in enumerate(B):
        for j, col in enumerate(columns):
            assert sum(x * y for x, y in zip(row, col)) == D * C_rhs[i][j]
    g = D
    for row in Q:
        for x in row:
            g = math.gcd(g, x)
    assert g == 1
    return {'BQ_equals_DC': True, 'positive_minimal_common_denominator': True,
            'arithmetic': 'Python arbitrary-precision integers',
            'nonsingularity': 'must be supplied by a separate exact certificate'}
