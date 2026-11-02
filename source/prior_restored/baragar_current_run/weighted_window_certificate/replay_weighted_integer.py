#!/usr/bin/env python3
"""Independent integer replay of exponential-window min(A,B/t), with direct integer-power rho checks.

Extends the frozen 164 integer checker, re-derives every harmonic floor/ceiling
prefix from every original node, and validates all splitting thresholds.

No imports from either producing program and no Arb requirement. Root bounds
are checked by integer powers; log bounds by an explicit positive atanh series;
pi bounds by Machin's identity and alternating-series remainder bounds. Every
old and new node is read, every block/prefix sum checked, and all 64 complete
allowed-prime recurrences are recomputed from the infinite-tail constant 1.
The direct analytic node enclosures remain an explicitly identified Arb input.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import struct
import sys

sys.set_int_max_str_digits(0)
BITS = 160
S = 1 << BITS
BASE = Path(__file__).resolve().parents[1]


def cd(a, b):
    assert a >= 0 and b > 0
    return (a + b - 1) // b


def mul(a, b):
    return cd(a * b, S)


def f(x):
    return F(int(x['numerator']), int(x['denominator']))


def upper(x):
    return cd(x.numerator * S, x.denominator)


def lower(x):
    return x.numerator * S // x.denominator


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        while b := stream.read(1024 * 1024):
            h.update(b)
    return h.hexdigest()


def atanh_log_bounds(z, count=140):
    assert 0 <= z <= F(1, 3)
    term, total, zz = z, F(0), z * z
    for j in range(count):
        total += 2 * term / (2 * j + 1)
        term *= zz
    remainder = 2 * term / ((2 * count + 1) * (1 - zz))
    return total, total + remainder


def atan_bounds(x, count):
    term, total = x, F(0)
    xx = x * x
    for j in range(count):
        total += term / (2 * j + 1)
        term *= -xx
    next_term = term / (2 * count + 1)
    return min(total, total + next_term), max(total, total + next_term)


def rootfloor(n, k):
    x = 1 << cd(n.bit_length(), k)
    while True:
        y = ((k - 1) * x + n // (x ** (k - 1))) // k
        if y >= x:
            break
        x = y
    while (x + 1) ** k <= n:
        x += 1
    while x ** k > n:
        x -= 1
    assert x ** k <= n < (x + 1) ** k
    return x


def exact_parameter_checks(p):
    assert p['epsilon'] == '4/45' and p['sigma'] == '1/18'
    ps = p['active_primes']
    actual = [q for q in range(2, 431)
              if all(q % d for d in range(2, isqrt(q) + 1))]
    assert [x['q'] for x in ps] == actual
    a5l, a5u = atan_bounds(F(1, 5), 110)
    a239l, a239u = atan_bounds(F(1, 239), 32)
    piL, piU = 16 * a5l - 4 * a239u, 16 * a5u - 4 * a239l
    assert f(p['pi_lower']) <= piL <= piU <= f(p['pi_upper'])
    ln2L, ln2U = atanh_log_bounds(F(1, 3))
    candidate_checks = 0
    for row in ps:
        q = row['q']
        for key, num, den in [('w', 4, 45), ('r', 1, 18)]:
            lo, hi = f(row[key + '_lower']), f(row[key + '_upper'])
            assert 0 < lo <= hi < 1
            assert lo.numerator ** den * q ** num <= lo.denominator ** den
            assert hi.numerator ** den * q ** num >= hi.denominator ** den
        k = q.bit_length() - 1
        z = F(q - (1 << k), q + (1 << k))
        zL, zU = atanh_log_bounds(z)
        lo, hi = k * ln2L + zL, k * ln2U + zU
        assert f(row['logq_lower']) <= lo <= hi <= f(row['logq_upper'])
        w, r = f(row['w_upper']), f(row['r_upper'])
        # Deliberately coarser 160-bit endpoints admit an exact rational check
        # against the outward root endpoints without relying on Arb K/L balls.
        K = F(upper(f(row['K_upper'])), S)
        L = F(upper(f(row['L_upper'])), S)
        rp, s0, s1 = F(1), F(1), F(0)
        E = row['first_discarded_exponent']
        for e in range(1, E):
            rp *= r
            s0 += rp
            s1 += e * rp
            assert w ** e * s0 <= K and w ** e * s1 <= L
            candidate_checks += 1
        assert w ** E / (1 - r) < 1
    # Independent rational roots at the first omitted prime, no Arb input.
    R = 1 << 512
    wu = F(rootfloor(R ** 45 // 431 ** 4, 45) + 1, R)
    ru = F(rootfloor(R ** 18 // 431, 18) + 1, R)
    assert wu ** 45 * 431 ** 4 >= 1 and ru ** 18 * 431 >= 1
    assert wu * (1 + ru) < 1 and p['first_omitted_prime'] == 431
    return {'prime_count': len(ps), 'polynomial_candidates_checked': candidate_checks,
        'all_root_intervals_by_integer_powers': True,
        'all_log_intervals_by_positive_series': True,
        'pi_interval_by_Machin_series': True,
        'K_L_exact_checks_after_160_bit_upward_rounding': True,
        'first_omitted_prime_independently_verified': 431}


def read_nodes(path, columns, N, blocks, blockkey, rawhash, prefixes=None, harmonic=None):
    codec = struct.Struct('<' + 'Q' * columns)
    totals, bs, bi, pi = [0] * columns, [0] * columns, 0, 0
    h = hashlib.sha256()
    hlow, hupp, hi = [0] * columns, [0] * columns, 0
    hrecords, hbits = (harmonic['records'], harmonic['metadata']['harmonic_scale_bits']) if harmonic else (None, None)
    with gzip.open(path, 'rb') as stream:
        for node in range(N):
            b = stream.read(codec.size)
            assert len(b) == codec.size
            h.update(b)
            for j, v in enumerate(codec.unpack(b)):
                totals[j] += v
                bs[j] += v
                if harmonic and node:
                    numerator = v * (1 << hbits)
                    hlow[j] += numerator // node
                    hupp[j] += (numerator + node - 1) // node
            if harmonic and node + 1 == hrecords[hi]['panels']:
                hp = hrecords[hi]
                assert totals == list(map(int, hp['node_sums']))
                assert hlow == list(map(int, hp['harmonic_lower_scaled']))
                assert hupp == list(map(int, hp['harmonic_upper_scaled']))
                hi += 1
                if hi == len(hrecords):
                    assert node + 1 == N
            if node == blocks[bi]['last_node']:
                assert bs == blocks[bi][blockkey]
                bs = [0] * columns
                bi += 1
            if prefixes and node + 1 == prefixes[pi]['panels']:
                assert totals == prefixes[pi]['node_sums']
                pi += 1
                if pi == len(prefixes):
                    assert node + 1 == N
        assert stream.read(1) == b'' and bi == len(blocks)
    assert h.hexdigest() == rawhash
    if prefixes:
        assert pi == len(prefixes)
    if harmonic:
        assert hi == len(hrecords)
    return totals


def integrals(d, path, refined):
    p, info, old = d['parameters'], d['integration'], d['legacy_rank_certificate']
    ps = p['active_primes']
    n, start, count = len(ps), 7, 64
    rankpath, suppath = BASE / old['rank_json'], BASE / old['support_json']
    assert sha(rankpath) == old['rank_json_sha256'] == 'c829047ffefa4a60c8b8c988dc07abb143a93985f3a33fe23852b2fdf4605871'
    assert sha(suppath) == old['support_json_sha256'] == 'b93e77f8e99d2ce4c930b6a3d02d068d7ee22ece2bd88e2a6d0a138a4477c4ac'
    rank, sup = json.loads(rankpath.read_text()), json.loads(suppath.read_text())
    ri, si = rank['integration'], sup['integration']
    rp, sp = rankpath.parent / ri['node_file'], suppath.parent / si['node_file']
    assert sha(rp) == ri['node_gzip_sha256'] == old['rank_node_gzip_sha256']
    assert sha(sp) == si['gzip_sha256'] == old['support_node_gzip_sha256']
    rs = read_nodes(rp, n + 1, ri['panels'], ri['blocks'], 'tail_sums', ri['node_raw_sha256'])
    ss = read_nodes(sp, count, si['panels'], si['blocks'], 'subset_sums', si['raw_sha256'])
    assert rs == [int(t['node_sum_integer']) for t in ri['tails']]
    assert ss == [int(t['node_sum_integer']) for t in si['subset_integrals']]
    newpath = path.parent / info['node_file']
    assert sha(newpath) == info['gzip_sha256']
    assert info['columns'] == n + 1 + count and info['tail_columns'] == n + 1
    assert info['panels'] == info['density'] * max(info['lambdas'])
    read_nodes(newpath, info['columns'], info['panels'], info['blocks'],
               'column_sums', info['raw_sha256'], info['prefix_integrals'], refined)
    K = [upper(f(t['K_upper'])) for t in ps]
    L = [upper(f(t['L_upper'])) for t in ps]
    logs = [upper(f(t['logq_upper'])) for t in ps]
    KK = [S] * (n + 1)
    for j in range(n - 1, -1, -1):
        KK[j] = mul(K[j], KK[j + 1])
    caps = KK.copy()
    for mask in range(count):
        cap = KK[start]
        for j in range(1, start):
            if mask & (1 << (j - 1)):
                cap = mul(cap, K[j])
        caps.append(cap)
    def corrections(factors):
        fs = [S] * (n + 1)
        for j in range(n - 1, -1, -1):
            fs[j] = mul(fs[j + 1], S + factors[j])
        for mask in range(count):
            cor = fs[start]
            for j in range(1, start):
                if mask & (1 << (j - 1)):
                    cor = mul(cor, S + factors[j])
            fs.append(cor)
        return fs
    N, scale = ri['panels'], 1 << ri['node_scale_bits']
    assert N == si['panels'] and ri['node_scale_bits'] == si['node_scale_bits']
    delta = cd(upper(f(p['pi_upper'])), 2 * N)
    fr = corrections([mul(x, delta) for x in L])
    rank_sums = rs + ss
    W = [min(caps[j], cd(v * fr[j], N * scale)) for j, v in enumerate(rank_sums)]
    density, scale = info['density'], 1 << info['node_scale_bits']
    fm = corrections([cd(mul(a, b), 2 * density) for a, b in zip(L, logs)])
    I = {}
    for pref in info['prefix_integrals']:
        lam = pref['lambda']
        I[lam] = [min(lam * caps[j], cd(v * fm[j], density * scale))
                  for j, v in enumerate(pref['node_sums'])]
    refined['correction_scaled'] = fm
    return K, W, I, {'old_rank_values': N * (n + 1), 'old_support_values': N * count,
        'all_harmonic_prefix_floors_and_ceilings_recomputed': True,
        'new_values': info['panels'] * info['columns'], 'all_blocks_and_prefixes_checked': True,
        'all_three_gzip_and_raw_hashes_checked': True}


def run_recursion(d, K, W, I, refined):
    p = d['parameters']
    ps, start = p['active_primes'], 7
    q, n = [t['q'] for t in ps], len(ps)
    w = [upper(f(t['w_upper'])) for t in ps]
    r = [upper(f(t['r_upper'])) for t in ps]
    logs = [upper(f(t['logq_upper'])) for t in ps]
    piL = lower(f(p['pi_lower']))
    assert piL > 3 * S
    meta, pref = refined['metadata'], refined['records']
    density, nodeScale = d['integration']['density'], 1 << d['integration']['node_scale_bits']
    Hscale = 1 << meta['harmonic_scale_bits']
    prefixes = {x['panels']: x for x in pref}
    thresholds = {(x['q'], x['a']): x for x in meta['thresholds']}
    universal = {x['q']: x for x in meta['universal_exponent_bounds']}

    @lru_cache(None)
    def refined_bound(j, a, column, lam):
        threshold = thresholds[q[j], a]
        N = lam * density
        k = min(N, threshold['split_panels'])
        end, begin = prefixes[N], prefixes[k]
        Sn, Sk = int(end['node_sums'][column]), int(begin['node_sums'][column])
        harmonic = (int(end['harmonic_upper_scaled'][column])
                    - int(begin['harmonic_lower_scaled'][column])) if k < N else 0
        assert harmonic >= 0
        Lupper = upper(f(threshold['L_upper']))
        Bupper = upper(f(threshold['B_upper']))
        aterm = cd(Lupper * Sk * S, piL * density * nodeScale)
        bterm = cd(Bupper * harmonic * S, piL * Hscale * nodeScale)
        cterm = cd(2 * Sn * S, lam * density * nodeScale)
        return mul(refined['correction_scaled'][column], aterm + bterm + cterm)

    @lru_cache(None)
    def coeff(j, a, lam):
        return cd(a * logs[j] * S, piL) + cd(2 * S, lam)

    @lru_cache(None)
    def window(j, a, index, column):
        kap = 1
        while q[index] ** kap < q[j] ** a:
            kap += 1
        best, method, lamchoice = kap * W[column], 'rank Sperner', None
        for lam, values in I.items():
            candidate = mul(coeff(j, a, lam), values[column])
            if candidate < best:
                best, method, lamchoice = candidate, 'Selberg Mellin', lam
            candidate = refined_bound(j, a, column, lam)
            if candidate < best:
                best, method, lamchoice = candidate, 'weighted exponential-window Selberg min(A,B/t)', lam
        return best, method, kap, lamchoice

    allC, oneC, records = [], [], []
    for mask in range(64):
        seq = [0] + [j for j in range(1, start) if mask & (1 << (j - 1))] + list(range(start, n))
        tails = [S] * (len(seq) + 1)
        columns = []
        for index in seq:
            if index >= start:
                columns.append(index)
            elif index > 0:
                suffix = sum(1 << (h - 1) for h in seq if index <= h < start)
                columns.append(n + 1 + suffix)
            else:
                columns.append(-1)
        for t in range(len(seq) - 1, -1, -1):
            tails[t] = mul(K[seq[t]], tails[t + 1])
        C, one, details = S, None, []
        for t in range(len(seq) - 1, -1, -1):
            j, nextC = seq[t], C
            assert C >= S and tails[t + 1] <= upper(f(universal[q[j]]['Euler_tail_cap']))
            wa, ra, a, best, besta, chosen = S, S, 0, C, 0, None
            while True:
                a += 1
                wa, ra = mul(wa, w[j]), mul(ra, r[j])
                stop = mul(wa, tails[t + 1] + C)
                if stop <= C:
                    if j == 0 and a == 1:
                        one = mul(wa, tails[t + 1] + mul(ra, C))
                    break
                J, E, method = tails[t + 1], S, ['Euler', None, None, None]
                for u in range(t + 1, len(seq)):
                    val, met, kap, lam = window(j, a, seq[u], columns[u])
                    cand = mul(E, val)
                    if cand < J:
                        J, method = cand, [met, q[seq[u]], kap, lam]
                    E = mul(E, K[seq[u]])
                value = mul(wa, J + mul(ra, C))
                if j == 0 and a == 1:
                    one = value
                if value > best:
                    best, besta, chosen = value, a, method
            assert a <= universal[q[j]]['first_omitted_exponent']
            C = best
            details.append({'q': q[j], 'C_next_scaled': str(nextC), 'C_scaled': str(C),
                'first_omitted_exponent': a, 'tail_test_scaled': str(stop),
                'maximizing_exponent': besta, 'maximizing_method': chosen})
        assert one is not None
        allC.append(C)
        oneC.append(one)
        records.append({'mask': mask, 'C_all_scaled': str(C), 'C_one_scaled': str(one), 'rows': details})
    joint = [allC[m] + oneC[63 ^ m] for m in range(64)]
    worst, top = max(enumerate(joint), key=lambda t: t[1])
    return cd(top, S), {'all_64_full_sequences_recomputed': True,
        'joint_upper_scaled': str(top), 'joint_integer_upper': cd(top, S),
        'worst_high_v2_mask': worst, 'assignment_upper_scaled': list(map(str, joint)),
        'support_records': records, 'cached_window_count': window.cache_info().currsize}


def endpoint(C, power):
    P = 10 ** power
    T = rootfloor((P + 1) ** 5, 44) + 1
    if (T - 1) ** 44 == (P + 1) ** 5:
        T -= 1
    lhs = (3*C)**6660 * (2*181265+3)**900 * (P+1)**673
    rhs = (16 * P) ** 900
    return {'power': power, 'T': str(T), 'main_exact': lhs < rhs,
        'root_exact': (T - 1) ** 44 < (P + 1) ** 5 <= T ** 44,
        'safe8': 36864 * 2541933 ** 4 * 2541935 ** 2 * T ** 6 < P,
        'primitive_cage': 256 * 2541933 ** 4 * T ** 4 < P}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrals', required=True)
    ap.add_argument('--recursion', required=True)
    ap.add_argument('--prefixes', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    ip, rp = Path(args.integrals), Path(args.recursion)
    d, first = json.loads(ip.read_text()), json.loads(rp.read_text())
    assert first['input_sha256'] == sha(ip)
    pars = exact_parameter_checks(d['parameters'])
    prefpath = Path(args.prefixes)
    meta = json.loads(prefpath.read_text())
    assert first['prefix_metadata_sha256'] == sha(prefpath) and meta['input_sha256'] == sha(ip)
    packed = prefpath.parent / meta['prefix_file']
    assert sha(packed) == meta['prefix_gzip_sha256']
    with gzip.open(packed, 'rb') as stream:
        raw = stream.read()
    assert hashlib.sha256(raw).hexdigest() == meta['prefix_raw_sha256']
    records = json.loads(raw)
    qp = {t['q']: t for t in d['parameters']['active_primes']}
    universal = {t['q']: t for t in meta['universal_exponent_bounds']}
    thresholds = {(t['q'], t['a']): t for t in meta['thresholds']}
    needed = {x['panels'] for x in d['integration']['prefix_integrals']}
    total = 0
    for q, row in qp.items():
        u = universal[q]
        A, cap = u['first_omitted_exponent'], f(u['Euler_tail_cap'])
        assert f(row['w_upper']) ** A * (cap + 1) <= 1
        for a in range(1, A):
            t = thresholds[q, a]
            rl, ru = f(t['rho_lower']), f(t['rho_upper'])
            assert 0 < rl <= ru < 1
            # Independent of the producer's iterative rounding: exact powers.
            assert rl.numerator**18*q**a <= rl.denominator**18
            assert ru.numerator**18*q**a >= ru.denominator**18
            L, B = 18*(1-rl), 1+ru
            cut = B*meta['density']/L
            split = cut.numerator//cut.denominator+1
            assert f(t['L_upper']) == L and f(t['B_upper']) == B
            assert t['split_panels'] == split
            needed.add(min(split, meta['panels']))
            total += 1
    assert total == len(thresholds) == len(meta['thresholds'])
    assert sorted(needed) == [t['panels'] for t in records]
    refined = {'metadata': meta, 'records': records}
    K, W, I, reading = integrals(d, ip, refined)
    C, rec = run_recursion(d, K, W, I, refined)
    assert C == first['recursion']['joint_integer_upper']
    power = first['endpoint']['power']
    ep, low = endpoint(C, power), endpoint(C, power - 1)
    okay = all(ep[k] for k in ['main_exact', 'root_exact', 'safe8', 'primitive_cage'])
    result = {'status': 'PASSED' if okay else 'ENDPOINT FAILED',
        'integer_fraction_only': True, 'fixed_bits': BITS, 'parameter_checks': pars,
        'node_readback': reading, 'integrals_sha256': sha(ip), 'first_recursion_sha256': sha(rp),
        'script_sha256': sha(Path(__file__)), 'same_integer_constant': True,
        'prefix_metadata_sha256': sha(prefpath), 'all_thresholds_verified': True,
        'all_exponential_weight_intervals_by_direct_integer_powers': True,
        'recursion': rec, 'endpoint': ep, 'lower_endpoint_diagnostic': low}
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'recursion'} |
        {'C': C, 'worst_mask': rec['worst_high_v2_mask']}, indent=2), flush=True)
    if not okay:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
