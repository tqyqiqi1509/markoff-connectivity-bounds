#!/usr/bin/env python3
"""Root closeout for the already reviewed prime-d general-rank proof chain.

The heavy producer and independent modular replay must have finished. This
binds their bytes, reads every journal cell and independently checks coverage.
It is not a replacement for that heavy replay or the rank-to-graph proof.
"""
import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
VERIFIER_SHA = 'b8d8a0f34c95df4f8dc1e2a0b87c7ee494050542769f12f41351d53fbc8e2f5e'
SERIAL_VERIFIER_SHA = 'b7f0ec03ab61b29e2976172722b8af022d2f43c332eed92db7989cfb62c5c424'
IMPLEMENTATION_ID = 'proven_proth_parallel_nmod_v1'
PRODUCER_SHA = 'e66db1d3ae87791a0eb5a55fae394bb0c33e69e16020ecf0f59755850bc555bc'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def accept(d, directory, replay_directory, commit=False):
    assert d > 2 and all(d % q for q in range(2, math.isqrt(d)+1))
    directory = directory.resolve()
    replay_directory = replay_directory.resolve()
    assert directory != replay_directory
    cp = directory / f'd{d}_universal_certificate.json'
    vp = replay_directory / f'd{d}_independent_universal_verification.json'
    hp = replay_directory / f'd{d}_modular_det_header.json'
    jp = replay_directory / f'd{d}_modular_det_journal.jsonl'
    c, v, h = [json.loads(p.read_text()) for p in (cp, vp, hp)]
    source = directory / c['source_file']
    code = Path(__file__).resolve().parent / 'verify_universal_modular_parallel.py'
    serial_code = directory / 'verify_universal_modular.py'
    assert sha(serial_code) == SERIAL_VERIFIER_SHA
    assert v['implementation_id'] == h['implementation_id'] == IMPLEMENTATION_ID
    assert v['serial_predecessor_sha256'] == h['serial_predecessor_sha256'] == SERIAL_VERIFIER_SHA
    assert type(v['parallel_workers']) is int and 1 <= v['parallel_workers'] <= 8
    producer = directory / 'universal_codimension_certificate.py'
    assert sha(code) == v['verifier_sha256'] == h['code_sha256'] == VERIFIER_SHA
    assert sha(producer) == c['producer_code_sha256'] == PRODUCER_SHA
    assert sha(source) == c['source_sha256'] == v['source_sha256'] == h['source_sha256']
    assert sha(cp) == v['certificate_sha256'] == h['certificate_sha256']
    assert sha(jp) == v['journal_sha256']
    assert c['d'] == v['d'] == h['d'] == d and c['N'] == v['N'] == 4*d
    assert c['bound'] == v['bound'] == 8*d
    assert c['p_conditions'] == {'strict_lower_bound': 8*d, 'excluded_residues_modulus': 2*d, 'excluded_residues': [1, 2*d-1]}
    assert c['status'] == 'CERTIFIED_ALL_ELIGIBLE_PRIMES'
    assert v['status'] == 'UNIVERSAL_CERTIFICATE_INDEPENDENTLY_VERIFIED'
    for key in ['all_raw_rows_reconstructed_and_normalizations_checked', 'primitive_lower_kernel_verified',
                'lower_maximal_minor_gcd_smooth', 'top_integer_system_verified_after_common_denominator',
                'residual_rough_gcd_one', 'integer_determinants_verified_by_explicit_modular_bound']:
        assert v[key] is True
    H = {k: int(x) for k, x in v['Hadamard_bounds'].items()}
    assert set(H) == {'lower', 'top'}
    D = {'lower': int(c['lower_certificate']['first_square_minor_det']), 'top': int(c['top_block_determinant'])}
    for name in H:
        assert 0 < abs(D[name]) <= H[name]
        assert H[name].bit_length() == h[name+'_H_bits']
    count = {k: 0 for k in H}
    covered = {k: 1 for k in H}
    product = 1
    last_q = 0
    first_done = {}
    lines = 0
    with jp.open() as f:
        for lines, line in enumerate(f, 1):
            assert line.endswith('\n')
            item = json.loads(line)
            q, k = item['q'], item['k']
            assert type(q) is int and type(k) is int and q > last_q
            assert 0 < k < 2**32 and k % 2 == 1 and q == k*2**32+1 and q < 2**64
            assert item['proth_exponent'] == 32 and item['witness'] == 3
            assert pow(3, (q-1)//2, q) == q-1
            last_q = q
            for name in H:
                value = item[name+'_det_mod_q']
                if product <= 2*H[name]:
                    assert type(value) is int and value == D[name] % q
                    count[name] += 1
                    covered[name] *= q
                else:
                    assert value is None
            product *= q
            for name in H:
                if name not in first_done and product > 2*H[name]:
                    first_done[name] = product
    assert lines == v['proven_moduli_count'] and count == v['checks_per_matrix']
    assert first_done == {k: int(x) for k, x in v['first_sufficient_modulus_products'].items()}
    assert all(covered[k] > 2*H[k] for k in H)
    source_state_path = directory/'STATE.json'
    source_state = json.loads(source_state_path.read_text())
    assert source_state['d'] == d
    assert (source_state['stage'] == 'STOPPED_AFTER_PHASE' and source_state.get('next_phase') == 'BOUNDED_MODULAR_REPLAY') or source_state['stage'] == 'COMPLETE_PENDING_ROOT_ACCEPTANCE'
    state_path = replay_directory/'PARALLEL_STATE.json'
    state = json.loads(state_path.read_text())
    assert state['d'] == d and state['stage'] == 'PARALLEL_COMPLETE_PENDING_ROOT_ACCEPTANCE'
    assert state['implementation_id'] == IMPLEMENTATION_ID and state['verifier_sha256'] == VERIFIER_SHA
    assert state['source_sha256'] == v['source_sha256'] and state['certificate_sha256'] == sha(cp)
    assert state['verification_sha256'] == sha(vp) and state['journal_sha256'] == sha(jp)
    assert state['moduli'] == lines and state['workers'] == v['parallel_workers']
    paths = [cp, vp, hp, jp, source, code, serial_code, producer, source_state_path, state_path]
    out = {
        'status': 'ACCEPTED_GENERAL_PRIME_D_CRITERION', 'accepted_utc': datetime.now(timezone.utc).isoformat(),
        'd': d, 'N': 4*d, 'condition': f'p prime, p>{8*d}, p mod {2*d} not in {{1,{2*d-1}}}',
        'conclusion': 'Pure-three-Vieta connectivity of all nonorigin points of x^2+y^2+z^2=xyz over F_p.',
        'scope': 'Reviewed general prime-d rank-to-graph theorem and pure-Vieta bridge are inherited dependencies; this is not an all-d saturation theorem.',
        'root_review': 'Explicitly separate parallel execution variant with pinned SHA. Source/algebra audit unchanged from the reviewed serial verifier. Fresh full modular replay completed under root ownership; every journal cell, proven-prime witness, per-matrix bound and final source/certificate/code hash rechecked here.',
        'parallel_implementation_id': IMPLEMENTATION_ID, 'parallel_workers': v['parallel_workers'],
        'serial_predecessor_sha256': SERIAL_VERIFIER_SHA,
        'modular_primes': lines, 'modular_determinant_checks': count,
        'actual_covered_modulus_product_bits': {k: covered[k].bit_length() for k in H},
        'hadamard_bits': {k: H[k].bit_length() for k in H},
        'all_journal_cells_present': True,
        'arithmetic_independence': 'Integer fmpz_mat and modular nmod_mat paths share FLINT; source/certificate logic and explicit bound checks are separate.',
        'source_reconstruction_scope': 'Every stored raw coefficient and normalization checked against the matrix; no claim every C++ coefficient was freshly generated by a second independent implementation.',
        'bindings': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'root_closeout_code_sha256': sha(Path(__file__).resolve()),
    }
    if not commit:
        print(json.dumps({'status': 'ALL_PARALLEL_CLOSEOUT_GATES_PASS_NO_ACCEPTANCE_WRITE', 'd': d, 'journal_lines': lines, 'checks': count, 'verifier_sha256': VERIFIER_SHA}), flush=True)
        return
    target = ROOT / f'control/D{d}_ROOT_ACCEPTANCE.json'
    assert not target.exists(), 'Acceptance is immutable: use an explicit superseding file if needed.'
    target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
    status_path = ROOT / 'control/CURRENT_STATUS.json'
    current = json.loads(status_path.read_text())
    current['accepted_new_general_criteria'] = sorted(set(current['accepted_new_general_criteria']+[d]))
    current['pending_general_parameters'] = [x for x in current['pending_general_parameters'] if x != d]
    sid = current['current_root_sessions'].pop(f'uniform{d}', None)
    if sid is not None:
        current.setdefault('completed_root_sessions', {})[f'uniform{d}'] = sid
    parallel_sid = current['current_root_sessions'].pop(f'uniform{d}_parallel', None)
    if parallel_sid is not None:
        current.setdefault('completed_root_sessions', {})[f'uniform{d}_parallel'] = parallel_sid
    current['updated_utc'] = current['last_updated_utc'] = out['accepted_utc']
    status_path.write_text(json.dumps(current, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'d': d, 'status': out['status'], 'acceptance_sha256': sha(target), 'journal_lines': lines, 'checks': count}), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('d', type=int)
    ap.add_argument('directory', type=Path)
    ap.add_argument('replay_directory', type=Path)
    ap.add_argument('--commit', action='store_true', help='Root-only explicit final acceptance after full review; otherwise read-only gates.')
    args = ap.parse_args()
    accept(args.d, args.directory, args.replay_directory, args.commit)
