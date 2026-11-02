#!/usr/bin/env python3
"""Exact assembly of accepted rank filters with the analytic support split.

No new determinant or spectral-node evaluation is claimed. Every reused
certificate is bound by SHA256; all new combinations use Python integers.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
from math import gcd, lcm, prod
from datetime import datetime, timezone
import hashlib
import json
import sys

sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACK = ROOT / 'agents/attachment_audit/extracted'
BINDINGS = {}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p, expected=None):
    p = ROOT / p if not isinstance(p, Path) else p
    digest = sha(p)
    if expected is not None:
        assert digest == expected, (str(p), digest, expected)
    BINDINGS[str(p.relative_to(ROOT))] = digest
    return json.loads(p.read_text())


def frac(j):
    return F(int(j['numerator']), int(j['denominator']))


def pack(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def factors(n):
    ans = {}
    q = 2
    while q*q <= n:
        while n % q == 0:
            ans[q] = ans.get(q, 0) + 1
            n //= q
        q += 1
    if n > 1:
        ans[n] = ans.get(n, 0) + 1
    return ans


def main():
    ds = [22,33,15,35,39,34,38,46,58]
    gate_rows = []
    for d in ds:
        if d == 58:
            acc = read('control/D58_ROOT_ACCEPTANCE.json',
                       'b4a20cbc0fedf5f6a5f67c5aa64721dc1022f3b2110c99d9331701f9af27323a')
            assert acc['status'] == 'ACCEPTED_GENERAL_D58_CRITERION'
            p = ROOT / 'agents/attachment_audit/d58/composite_d58_certificate.json'
            cert = read(p, acc['input_bindings'][str(p.relative_to(ROOT))])
            review = read(ROOT / 'agents/attachment_audit/d58/composite_d58_independent_verification.json',
                          acc['input_bindings']['agents/attachment_audit/d58/composite_d58_independent_verification.json'])
        elif d == 46:
            pbase = ROOT / 'prior_restored/latest/finite46_recomputed'
            acc = read(pbase / 'd46_recomputed_acceptance.json')
            assert acc['status'] == 'D46_FRESH_GENERAL_CERTIFICATE_AND_BOUNDED_MODULAR_REPLAY_COMPLETED'
            cert = read(pbase / 'composite_d46_certificate.json', acc['core_files']['composite_d46_certificate.json'])
            review = read(pbase / 'composite_d46_independent_verification.json', acc['core_files']['composite_d46_independent_verification.json'])
        else:
            p, = (PACK / '06_computations').glob(f'*_composite_d{d}_certificate.json')
            cert = read(p)
            if d in [15,35]:
                base = ROOT / 'prior_restored/original_v14/markoff_project/agents/computational_mechanisms'
                fn = ('composite_d15_independent_bareiss_verification.jsonl' if d == 15
                      else 'composite_d35_minor_verification.jsonl')
                review = read(base / fn)
                assert review['arithmetic_verified'] is True
                assert review['matrix_sha256'] == cert['matrix_sha256']
            else:
                rp, = (PACK / '04_audits_history').glob(f'*_composite_d{d}_independent_verification.json')
                review = read(rp)
                assert review['certificate_sha256'] == sha(p)
        assert cert['status'] == 'certified_no_exception_p_above_bound'
        assert int(cert['rough_gcd']) == 1
        N = int(cert['nd'])
        assert N >= 3*d
        if d not in [15,35]:
            assert review['status'] == 'GENERAL_COMPOSITE_UNIFORM_CERTIFICATE_INDEPENDENTLY_VERIFIED'
            assert review['strict_p_bound'] == 2*N
            assert review['excluded_residues_modulus'] == 2*d
            assert review['excluded_residues'] == [1,2*d-1]
            assert int(review['rough_gcd_of_verified_integer_minors']) == 1
        gate_rows.append({'d':d,'N':N,'strict_p_bound':2*N,'modulus':2*d})
    assert max(r['strict_p_bound'] for r in gate_rows[:6]) == 468
    assert max(r['strict_p_bound'] for r in gate_rows) == 696

    moduli = [2*d for d in ds]
    cumulative = [moduli[0]]
    overlaps = []
    for m in moduli[1:]:
        overlaps.append(gcd(cumulative[-1],m))
        cumulative.append(lcm(cumulative[-1],m))
    assert all(g > 2 for g in overlaps)
    compatible = []
    for signs in product([-1,1],repeat=len(ds)):
        if all((signs[i]-signs[j]) % gcd(moduli[i],moduli[j]) == 0
               for i,j in combinations(range(len(ds)),2)):
            compatible.append(signs)
    assert compatible == [(-1,)*9,(1,)*9]
    assert cumulative[5] == 1021020
    assert cumulative[-1] == 12939386460
    prime_support = sorted(factors(cumulative[-1]))
    assert prime_support == [2,3,5,7,11,13,17,19,23,29]

    # Independently rebuild the prime-power data for the already accepted52 set.
    all_ds = [4,5,7,8,9,11,13,15,16,17,19,22,23,25,27,29,31,32,33,34,35,37,38,39,41,43,46,47,49,53,58,59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,131,137]
    qpow = factors(lcm(*(2*d for d in all_ds)))
    forced_high_core = prod(q**(qpow[q]-(q==2)) for q in prime_support)
    assert forced_high_core == 10760393780136000
    finite = read('control/FINITE_1E56_ROOT_ACCEPTANCE.json',
                  '4e62c14593d5aa980f0a3645913cd4c390a828fdbd90697763ce59dd46bf734d')
    assert finite['status'] == 'ACCEPTED_INTERNAL_PROOF_CHAIN'
    enum = read('agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json',
                finite['bindings']['agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json'])
    assert enum['sign_components'][0] == prime_support
    assert {int(q):a for q,a in enum['prime_exponents'].items()} == qpow

    half_acc = read('control/HALF_BOUND_ROOT_ACCEPTANCE.json')
    hname = 'agents/joint_analytic/HALF_MODULUS_ENDPOINT_CERTIFICATE.json'
    half = read(hname, half_acc['input_sha256'][hname])
    assert half['trace']['first_inactive_prime'] == 401
    # Recombine actual accepted local upper bounds; never divide two upper bounds.
    trace_rows = [r for r in half['trace']['rows'] if r['q'] >= 31]
    low_B = prod((F(r['local_upper_integer'],1<<80) for r in trace_rows),start=F(1))
    assert low_B < 662
    neg_rows = [r for r in half['negative_odd']['local_rows'] if r['q'] >= 31]
    low_ratio = prod((frac(r['local_upper']) for r in neg_rows),start=F(1))
    assert half['negative_odd']['prime_tail_from'] == 73
    assert low_ratio < F(26,25)

    trace_unchanged = []
    for r in half['trace']['rows']:
        q = r['q']
        if q not in prime_support:
            continue
        amin = qpow[q]-(q==2)
        arg = max(r['finite_values'],key=lambda v:v['upper_integer'])
        allowed_upper = max([1<<80]+[v['upper_integer'] for v in r['finite_values'] if v['e']>=amin])
        assert allowed_upper == r['local_upper_integer'] and arg['e'] >= amin
        trace_unchanged.append({'q':q,'minimum_exponent':amin,
                                'existing_Euler_maximizing_exponent':arg['e'],
                                'existing_cap_unchanged':True})

    old63 = read('prior_restored/latest/root/bad_branch_transfer/old63branches_and_prerequisites.json')
    old_hashes = {Path(x['path']).name:x['sha256'] for x in old63['sources']}
    oldbase = ROOT / 'prior_restored/baragar_current_run/weighted_window_certificate'
    ordinary = read(oldbase/'weighted_159.json',old_hashes['weighted_159.json'])
    replay = read(oldbase/'replay_weighted_159.json',old_hashes['replay_weighted_159.json'])
    assert ordinary['status'] == 'ARITHMETIC PASSED' and replay['status'] == 'PASSED'
    pr = next(r for r in ordinary['recursion']['big_rows'] if r['q'] == 31)
    prod_count = frac(pr['C_upper'])
    rc_values = []
    for record in replay['recursion']['support_records']:
        rr = next(r for r in record['rows'] if r['q'] == 31)
        rc_values.append(F(int(rr['C_scaled']),1<<replay['fixed_bits']))
    assert len(rc_values) == 64 and len(set(rc_values)) == 1
    assert prod_count < 135 and rc_values[0] < 135

    # Check that mere forced minimum exponents do not lower the accepted BE K.
    be = read('agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json',
              'ca42290f5f622d7e71ab7f68041d2467350e45a0a91a519c8f08144daf8f9992')
    unchanged = []
    for r in be['euler']['rows']:
        q = r['q']
        if q not in prime_support:
            continue
        amin = qpow[q]-(q==2)
        assert r['maximizer_of_enclosures'] >= amin
        allowed_upper = max([F(1)]+[F(v['local_value_upper']['upper'])
                                    for v in r['values'] if v['a']>=amin])
        assert allowed_upper == F(r['local_upper']['upper'])
        unchanged.append({'q':q,'minimum_exponent':amin,
                          'existing_Euler_maximizing_exponent':r['maximizer_of_enclosures'],
                          'existing_cap_unchanged':True})

    result = {
      'status':'EXACT_SUPPORT_ASSEMBLY_PASSED_PENDING_ROOT_REVIEW',
      'utc':datetime.now(timezone.utc).isoformat(),
      'conjecture':'OPEN; no new numerical p endpoint asserted',
      'rank_gates':gate_rows,
      'successive_moduli':cumulative,
      'overlap_gcd':overlaps,
      'compatible_sign_patterns':[list(s) for s in compatible],
      'six_filter_modulus':1021020,'six_filter_strict_p_gate':468,
      'nine_filter_modulus':12939386460,'nine_filter_strict_p_gate':696,
      'forced_high_half_core_from52':str(forced_high_core),
      'forced_prime_support':prime_support,
      'low_half_prime_divisors_at_least':31,
      'eliminated_impossible_analytic_support_masks':list(range(63)),
      'remaining_analytic_support_mask':63,
      'low_half_constants':{
        'B_trace_cap':662,'B_trace_actual_product_upper':pack(low_B),
        'count_cap':135,'count_producer_upper':pack(prod_count),
        'count_independent_replay_upper':pack(rc_values[0]),
        'negative_ratio_cap':pack(F(26,25)),
        'negative_ratio_actual_product_upper':pack(low_ratio),
        'trace_factor_prime_rows':[r['q'] for r in trace_rows],
        'negative_factor_prime_rows':[r['q'] for r in neg_rows]},
      'forced_exponents_alone_do_not_improve_current_BE_cap':unchanged,
      'forced_exponents_alone_do_not_improve_separated_joint_B_cap':trace_unchanged,
      'bindings':BINDINGS,
      'script_sha256':sha(Path(__file__)),
      'limits':[
        'No determinant, source-matrix or spectral node rerun in this assembly.',
        'The old 7e153 and 3e153 no-filter proofs remain separate and unchanged.',
        'Arithmetic mask elimination applies only after using the accepted rank-filter dichotomy.',
        'No claims that the remaining arithmetic progression contains no primes.'
      ]}
    (HERE/'SUPPORT_COUPLING_CERTIFICATE.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','nine_filter_modulus','forced_high_half_core_from52','low_half_constants']},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
