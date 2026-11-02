#!/usr/bin/env python3
"""Deterministic, restartable integer determinant with proved word primes.

Independent optional engine. Does not edit production matrices/certificates.
Usage: python resumable_crt_det.py MATRIX.json --checkpoint STATE.json
Choose the same minor as gcd_minors with --minor-index j, or explicit row
indices with --row-indices JSON. Resume is automatic and source SHA bound.
--max-primes permits a controlled partial run; no partial run is a certificate.
"""
from __future__ import annotations
import argparse, hashlib, json, math, pathlib, sys, time

sys.set_int_max_str_digits(0)
HERE = pathlib.Path(__file__).resolve().parent
COMPUTE = HERE.parent / 'computational_mechanisms'
sys.path.insert(0, str(COMPUTE / 'deps'))
from flint import fmpz_mat, nmod_mat

SCHEMA = 1
PROTH_R = 32
INITIAL_K = 2**29 - 1
TRIAL_PRIMES = (3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)


def atomic_json(path: pathlib.Path, data):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, separators=(',', ':')) + '\n')
    temporary.replace(path)


def verified_proth(p: int, k: int, r: int, witness: int) -> bool:
    # This is a proof, not a probable-prime test. See PROOF_AND_COST.md.
    return (k > 0 and k % 2 == 1 and k < 2**r and
            p == k * 2**r + 1 and
            pow(witness, (p - 1)//2, p) == p - 1)


def next_proth(k: int):
    """Search deterministically downwards; accepting only proven primes."""
    while k >= 2**28 + 1:
        p = k * 2**PROTH_R + 1
        if all(p % q for q in TRIAL_PRIMES):
            for witness in range(2, 34):
                euler = pow(witness, (p - 1)//2, p)
                if euler == p - 1:
                    assert verified_proth(p, k, PROTH_R, witness)
                    return p, k, witness, k-2
                if euler != 1:
                    # Euler's necessary condition proves this candidate
                    # composite; do not waste the remaining witness trials.
                    break
        k -= 2
    raise RuntimeError('Proth search interval exhausted; do not silently reuse a prime')


def matrix_input(path: pathlib.Path, minor_index=0, row_indices=None):
    payload = path.read_bytes()
    source = json.loads(payload)
    allrows = source['rows'] if isinstance(source, dict) else source
    n = len(allrows[0]) if allrows else 0
    if not all(len(row) == n for row in allrows):
        raise ValueError('Ragged matrix')
    if len(allrows) < n:
        raise ValueError('Not enough rows for a square maximal minor')
    if row_indices is None:
        chosen = list(range(n))
        if minor_index:
            if len(allrows) == n:
                raise ValueError('No replacement row')
            chosen[(minor_index-1) % n] = n + (minor_index-1) % (len(allrows)-n)
    else:
        chosen = list(row_indices)
    if len(chosen) != n or len(set(chosen)) != n or any(i < 0 or i >= len(allrows) for i in chosen):
        raise ValueError('Invalid row selection')
    rows = [[int(x) for x in allrows[i]] for i in chosen]
    return rows, {'source_file': str(path.resolve()), 'source_sha256': hashlib.sha256(payload).hexdigest(),
                  'shape': [len(allrows), n], 'row_indices': chosen}


def content_reduce(rows):
    """Extract exact row/column contents; their product is restored at the end."""
    a = [row[:] for row in rows]
    n = len(a)
    factor = 1
    row_contents = []
    for i in range(n):
        g = math.gcd(*a[i])
        if g == 0:
            return a, 0, {'zero_row': i, 'row_contents': row_contents, 'column_contents': []}
        row_contents.append(g)
        factor *= g
        a[i] = [x//g for x in a[i]]
    col_contents = []
    for j in range(n):
        g = math.gcd(*(a[i][j] for i in range(n)))
        if g == 0:
            return a, 0, {'zero_column': j, 'row_contents': row_contents, 'column_contents': col_contents}
        col_contents.append(g)
        factor *= g
        for i in range(n):
            a[i][j] //= g
    return a, factor, {'row_contents': row_contents, 'column_contents': col_contents}


def hadamard_data(rows):
    n = len(rows)
    row_square = math.prod(sum(x*x for x in row) for row in rows)
    column_square = math.prod(sum(rows[i][j]**2 for i in range(n)) for j in range(n))
    h2 = min(row_square, column_square)
    # B=2^bound_bits is an exact upper bound on |det|. No floating logarithm.
    bound_bits = (h2.bit_length()+1)//2
    return {'square_bit_length': h2.bit_length(),
            'square_sha256_hex': hashlib.sha256(hex(h2).encode()).hexdigest(),
            'bound_bits': bound_bits,
            'proth_count_upper_60bit': (bound_bits+2+59)//60}


def crt_step(residue, modulus, p, determinant_mod_p):
    correction = ((determinant_mod_p-residue) * pow(modulus % p, -1, p)) % p
    return residue + modulus*correction, modulus*p


def validate_state(state, binding, reduction, hadamard, recompute=False, A=None):
    assert state['schema'] == SCHEMA
    # Environment recovery can move an identical source file. The path is
    # diagnostic; bytes, dimensions, and ordered row indices bind the proof.
    assert ({k:v for k,v in state['binding'].items() if k != 'source_file'} ==
            {k:v for k,v in binding.items() if k != 'source_file'})
    assert state['content_reduction'] == reduction
    assert state['hadamard'] == hadamard
    residue, modulus = 0, 1
    previous_k = INITIAL_K + 2
    for record in state['records']:
        p, k, r, witness, dm = (record[x] for x in ('p','k','r','witness','det_mod_p'))
        assert k < previous_k and verified_proth(p, k, r, witness)
        assert 0 <= dm < p
        if recompute:
            assert A is not None and int(nmod_mat(A, p).det()) == dm
        residue, modulus = crt_step(residue, modulus, p, dm)
        previous_k = k
    assert str(residue) == state['residue'] and str(modulus) == state['modulus']
    assert state['next_k'] <= previous_k-2
    return residue, modulus


def run(path, checkpoint, minor_index=0, row_indices=None, max_primes=None,
        max_seconds=None, checkpoint_every=32, verify_resume=False):
    started = time.monotonic()
    rows, binding = matrix_input(path, minor_index, row_indices)
    if not rows:
        return {'status': 'certified_exact', 'determinant': '1', 'binding': binding}
    reduced, factor, reduction = content_reduce(rows)
    if not factor:
        result = {'schema': SCHEMA, 'status': 'certified_exact_zero', 'determinant': '0',
                  'binding': binding, 'content_reduction': reduction}
        atomic_json(checkpoint, result)
        return result
    hadamard = hadamard_data(reduced)
    A = fmpz_mat(reduced)
    target = 1 << (hadamard['bound_bits']+1)
    if checkpoint.exists():
        state = json.loads(checkpoint.read_text())
        residue, modulus = validate_state(state, binding, reduction, hadamard,
                                          recompute=verify_resume, A=A)
        assert int(state['content_factor']) == factor
    else:
        residue, modulus = 0, 1
        state = {'schema': SCHEMA, 'status': 'partial_not_a_certificate', 'binding': binding,
                 'content_reduction': reduction, 'content_factor': str(factor),
                 'hadamard': hadamard, 'records': [], 'next_k': INITIAL_K,
                 'residue': '0', 'modulus': '1', 'elapsed_sessions_seconds': 0.0}
    old_elapsed = state.get('elapsed_sessions_seconds', 0.0)
    count_this_run = 0
    while modulus <= target:
        if ((max_primes is not None and count_this_run >= max_primes) or
            (max_seconds is not None and time.monotonic()-started >= max_seconds)):
            break
        prime_started = time.monotonic()
        p, k, witness, next_k = next_proth(state['next_k'])
        modular = nmod_mat(A, p)  # Native FLINT conversion; no Python entry loop.
        detmod = int(modular.det())
        residue, modulus = crt_step(residue, modulus, p, detmod)
        state['records'].append({'p': p, 'k': k, 'r': PROTH_R, 'witness': witness,
                                 'det_mod_p': detmod,
                                 'seconds': time.monotonic()-prime_started})
        state['next_k'] = next_k
        count_this_run += 1
        if count_this_run % checkpoint_every == 0:
            state.update(residue=str(residue), modulus=str(modulus),
                         elapsed_sessions_seconds=old_elapsed+time.monotonic()-started)
            atomic_json(checkpoint, state)
            print(json.dumps({'status': 'partial_not_a_certificate',
                              'primes': len(state['records']), 'modulus_bits': modulus.bit_length(),
                              'target_bits': target.bit_length(),
                              'session_seconds': time.monotonic()-started}), flush=True)
    state.update(residue=str(residue), modulus=str(modulus),
                 elapsed_sessions_seconds=old_elapsed+time.monotonic()-started)
    if modulus > target:
        centered = residue if 2*residue <= modulus else residue-modulus
        assert abs(centered) <= (1 << hadamard['bound_bits'])
        state.update(status='certified_exact', reduced_determinant=str(centered),
                     determinant=str(factor*centered), rigorous_termination='modulus > 2 * exact_power_of_two_Hadamard_bound')
    else:
        state['status'] = 'partial_not_a_certificate'
    atomic_json(checkpoint, state)
    return state


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix', type=pathlib.Path)
    parser.add_argument('--checkpoint', type=pathlib.Path, required=True)
    parser.add_argument('--minor-index', type=int, default=0)
    parser.add_argument('--row-indices', type=json.loads)
    parser.add_argument('--max-primes', type=int)
    parser.add_argument('--max-seconds', type=float)
    parser.add_argument('--checkpoint-every', type=int, default=32)
    parser.add_argument('--verify-resume', action='store_true')
    args = parser.parse_args()
    if args.checkpoint_every < 1:
        raise ValueError('checkpoint-every must be positive')
    result = run(args.matrix, args.checkpoint, args.minor_index, args.row_indices,
                 args.max_primes, args.max_seconds, args.checkpoint_every, args.verify_resume)
    print(json.dumps({'status': result['status'],
                      'primes': len(result.get('records', [])),
                      'determinant_digits': len(result.get('determinant', '').lstrip('-')),
                      'checkpoint': str(args.checkpoint),
                      'elapsed_sessions_seconds': result.get('elapsed_sessions_seconds')}), flush=True)


if __name__ == '__main__':
    main()
