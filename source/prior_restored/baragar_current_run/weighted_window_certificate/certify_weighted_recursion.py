#!/usr/bin/env python3
"""Exact rational weighted-window Selberg/rank/Euler recurrences on all 64 supports.

Extends the frozen 164 recurrence with the independently derived nonnegative
error estimate |Mhat(t)| <= min(L,2/|t|)+2*pi/lambda and two-sided node prefixes.

This reads a full Mellin integral certificate. Analytic majorant existence and
the accepted geometric Markoff inputs are explicit external dependencies.
All arithmetic decisions here are rational/integer, not displayed decimals.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import gzip
import json
from pathlib import Path
import sys

from flint import arb, ctx

sys.set_int_max_str_digits(0)
S = 1 << 128


def ceilf(x):
    return -((-x.numerator) // x.denominator)


def up(x):
    return F(ceilf(x * S), S)


def pair(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def rat(x):
    return F(int(x['numerator']), int(x['denominator']))


def show(x):
    return (arb(x.numerator) / x.denominator).str(48)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rootceil(n, k):
    lo, hi = 0, 1
    while hi ** k < n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** k < n:
            lo = mid
        else:
            hi = mid
    assert lo ** k < n <= hi ** k
    return hi


def endpoint(C, power):
    P, KB = 10 ** power, 181265
    lhs = (3*C)**6660 * (2*KB+3)**900 * (P+1)**673
    rhs = (16 * P) ** 900
    T = rootceil((P + 1) ** 5, 44)
    safe = 36864 * 2541933 ** 4 * 2541935 ** 2 * T ** 6
    cage = 256 * 2541933 ** 4 * T ** 4
    ratio = (arb(3*C)**(arb(37)/5) * (arb(KB)+arb(3)/2) * arb(P+1)**(arb(673)/900) / (8*P))
    return {'power': power, 'C': C, 'K_B': KB, 'K_tau': 2541933, 'T': str(T),
        'main_left_sha256': hashlib.sha256(str(lhs).encode()).hexdigest(),
        'main_right_sha256': hashlib.sha256(str(rhs).encode()).hexdigest(),
        'checks': {'main_exact': lhs < rhs,
            'root_exact': (T - 1) ** 44 < (P + 1) ** 5 <= T ** 44,
            'safe8': safe < P, 'primitive_cage': cage < P},
        'main_ratio_display_only': ratio.str(48),
        'safe8_left': str(safe), 'primitive_cage_left': str(cage)}


def recurse(data, refined):
    p, info, old = data['parameters'], data['integration'], data['legacy_rank_certificate']
    rows, start, masks = p['active_primes'], 7, 64
    n = len(rows)
    q = [r['q'] for r in rows]
    assert n == 82 and q[1:start] == [3, 5, 7, 11, 13, 17]
    w = [rat(r['w_upper']) for r in rows]
    r = [rat(t['r_upper']) for t in rows]
    K = [rat(t['K_upper']) for t in rows]
    logs = [rat(t['logq_upper']) for t in rows]
    piL = rat(p['pi_lower'])
    assert 3 < piL < 4
    refmeta, pref = refined['metadata'], refined['records']
    Hscale = 1 << refmeta['harmonic_scale_bits']
    density, node_scale = info['density'], 1 << info['node_scale_bits']
    prefmap = {t['panels']: t for t in pref}
    thresholds = {(t['q'], t['a']): (t['split_panels'], rat(t['L_upper']), rat(t['B_upper']))
                  for t in refmeta['thresholds']}
    universal = {t['q']: t for t in refmeta['universal_exponent_bounds']}
    correction = [rat(t) for t in info['correction_upper']]

    @lru_cache(None)
    def refined_bound(j, a, column, lam):
        split, L, B = thresholds[q[j], a]
        N = lam * density
        k = min(split, N)
        end, beg = prefmap[N], prefmap[k]
        Sn, Sk = int(end['node_sums'][column]), int(beg['node_sums'][column])
        if k < N:
            harm = (int(end['harmonic_upper_scaled'][column])
                    - int(beg['harmonic_lower_scaled'][column]))
            assert harm >= 0
        else:
            harm = 0
        bound = (L / piL * F(Sk, density * node_scale)
                 + B / piL * F(harm, Hscale * node_scale)
                 + F(2 * Sn, lam * density * node_scale))
        return up(correction[column] * bound)
    waves = [rat(t) for t in old['rank_W_upper']] + [rat(t) for t in old['subset_W_upper']]
    mellin = {t['lambda']: [rat(v) for v in t['integral_upper']]
              for t in info['prefix_integrals']}
    KK = [F(1)] * (n + 1)
    for j in range(n - 1, -1, -1):
        KK[j] = up(K[j] * KK[j + 1])

    @lru_cache(None)
    def selberg_coefficient(j, a, lam):
        return up(a * logs[j] / piL + F(2, lam))

    @lru_cache(None)
    def interval(j, a, index, column):
        kap = 1
        qa = q[j] ** a
        while q[index] ** kap < qa:
            kap += 1
        assert q[index] ** kap >= qa and (kap == 1 or q[index] ** (kap - 1) < qa)
        best, method, lamchosen = up(kap * waves[column]), 'rank Sperner', None
        for lam, values in mellin.items():
            value = up(selberg_coefficient(j, a, lam) * values[column])
            if value < best:
                best, method, lamchosen = value, 'Selberg Mellin', lam
            value = refined_bound(j, a, column, lam)
            if value < best:
                best, method, lamchosen = value, 'weighted exponential-window Selberg min(A,B/t)', lam
        return best, method, kap, lamchosen

    def step(j, C, indices, columns, K0):
        assert C >= 1 and K0 <= rat(universal[q[j]]['Euler_tail_cap'])
        A = 1
        while w[j] ** A * (K0 + C) > C:
            A += 1
        stop = w[j] ** A * (K0 + C)
        assert stop <= C and A < 10000
        assert A <= universal[q[j]]['first_omitted_exponent']
        best, besta, one, tests = C, 0, None, []
        for a in range(1, A):
            J, method, cut, layers, lam = K0, 'Euler', None, None, None
            E = F(1)
            for index, column in zip(indices, columns):
                bound, met, kap, la = interval(j, a, index, column)
                trial = up(E * bound)
                if trial < J:
                    J, method, cut, layers, lam = trial, met, q[index], kap, la
                E = up(E * K[index])
            value = up(w[j] ** a * (J + r[j] ** a * C))
            tests.append({'a': a, 'J_upper': pair(J), 'method': method,
                'cut_prime': cut, 'layers': layers, 'lambda': lam, 'value_upper': pair(value)})
            if a == 1:
                one = value
            if value > best:
                best, besta = value, a
        # Fixed v2=1 is needed even if a=1 is already dominated by C.
        if one is None:
            assert A == 1
            one = up(w[j] * (K0 + r[j] * C))
        return best, one, {'q': q[j], 'C_next_upper': pair(C), 'K_tail_upper': pair(K0),
            'first_omitted_exponent': A, 'omitted_tail_test_upper': pair(stop),
            'maximizing_exponent': besta, 'C_upper': pair(best), 'candidates': tests}

    C, bigrows = F(1), []
    for j in range(n - 1, start - 1, -1):
        indices = list(range(j + 1, n))
        C, _, record = step(j, C, indices, indices, KK[j + 1])
        bigrows.append(record)
    Cbig = C
    constants, supportrows = [], []
    for mask in range(masks):
        included = [j for j in range(1, start) if mask & (1 << (j - 1))]
        def tail(j):
            rem = [h for h in included if h > j]
            indices = rem + list(range(start, n))
            columns = [n + 1 + sum(1 << (h - 1) for h in rem if h >= k) for k in rem]
            columns += list(range(start, n))
            K0 = F(1)
            for h in indices:
                K0 = up(K0 * K[h])
            return indices, columns, K0
        C, sr = Cbig, []
        for j in included[::-1]:
            C, _, record = step(j, C, *tail(j))
            sr.append(record)
        Codd = C
        Call, Cone, record = step(0, C, *tail(0))
        sr.append(record)
        constants.append({'mask': mask, 'allowed_small_primes': [q[j] for j in included],
            'C_odd_upper': pair(Codd), 'C_all_two_upper': pair(Call), 'C_v2_one_upper': pair(Cone),
            'C_all_two_display_only': show(Call), 'C_v2_one_display_only': show(Cone)})
        supportrows.append({'mask': mask, 'rows': sr})
    assignments = []
    for mask in range(masks):
        bound = rat(constants[mask]['C_all_two_upper']) + rat(constants[63 ^ mask]['C_v2_one_upper'])
        assignments.append({'high_v2_mask': mask, 'low_v2_mask': 63 ^ mask,
            'joint_upper': pair(bound), 'joint_display_only': show(bound)})
    worst = max(assignments, key=lambda d: rat(d['joint_upper']))
    joint = rat(worst['joint_upper'])
    return {'arithmetic': 'Exact Fraction, positive upward rounding to 2^-128',
        'big_constant_upper': pair(Cbig), 'big_rows': bigrows, 'support_rows': supportrows,
        'support_constants': constants, 'assignments': assignments, 'worst_assignment': worst,
        'joint_integer_upper': ceilf(joint), 'joint_display_only': show(joint),
        'interval_cache_entries': interval.cache_info().currsize}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--prefixes', required=True)
    ap.add_argument('--power', type=int, default=159)
    args = ap.parse_args()
    ctx.prec = 256
    source, out = Path(args.input), Path(args.output)
    data = json.loads(source.read_text())
    assert data['parameters']['epsilon'] == '4/45' and data['parameters']['sigma'] == '1/18'
    prefixpath = Path(args.prefixes)
    meta = json.loads(prefixpath.read_text())
    assert meta['input_sha256'] == sha(source)
    packed = prefixpath.parent / meta['prefix_file']
    assert sha(packed) == meta['prefix_gzip_sha256']
    with gzip.open(packed, 'rb') as stream:
        raw = stream.read()
    assert hashlib.sha256(raw).hexdigest() == meta['prefix_raw_sha256']
    result = recurse(data, {'metadata': meta, 'records': json.loads(raw)})
    C = result['joint_integer_upper']
    ep, lower = endpoint(C, args.power), endpoint(C, args.power - 1)
    status = 'ARITHMETIC PASSED' if all(ep['checks'].values()) else 'ENDPOINT FAILED'
    output = {'status': status, 'input_file': str(source), 'input_sha256': sha(source),
        'prefix_metadata_file': str(prefixpath), 'prefix_metadata_sha256': sha(prefixpath),
        'script_sha256': sha(Path(__file__)), 'recursion': result,
        'endpoint': ep, 'lower_endpoint_diagnostic': lower,
        'analytic_dependencies': ['positive exponential mixture of Selberg majorants; research/weighted_window_selberg_proof.md', 'geometric coefficient729/8; research/full_denominator_cutoff_proof.md',
            'previously audited weighted strong Sperner theorem',
            'previously audited Markoff cage, safe8, Wronskian and pure Vieta transfer']}
    out.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'status': status, 'C': C, 'joint': result['joint_display_only'],
        'worst_high_v2_mask': result['worst_assignment']['high_v2_mask'],
        'endpoint': ep, 'lower': lower}, indent=2), flush=True)
    if status != 'ARITHMETIC PASSED':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
