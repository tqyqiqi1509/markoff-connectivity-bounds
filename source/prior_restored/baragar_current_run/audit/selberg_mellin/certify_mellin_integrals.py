#!/usr/bin/env python3
"""Rigorous Mellin-frequency integrals for a Selberg interval majorant.

This program does not import the floating-point exploration or its recurrence.
It certifies every saved node by direct Arb/acb finite polynomial summation,
and integrates all 83 prime tails and all 64 small-prime supports. Old rank
nodes are accepted predecessor certificates, whose complete streams and all
block sums are rechecked. No output of this program alone audits Selberg's
analytic majorant theorem; that source/theorem check is a separate dependency.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import struct
import sys
import time

import flint
from flint import arb, acb, ctx

sys.set_int_max_str_digits(0)
BITS = 128
GRID = 1 << BITS
BASE = Path(__file__).resolve().parents[2]
EXPECTED_RANK = 'c829047ffefa4a60c8b8c988dc07abb143a93985f3a33fe23852b2fdf4605871'
EXPECTED_SUPPORT = 'b93e77f8e99d2ce4c930b6a3d02d068d7ee22ece2bd88e2a6d0a138a4477c4ac'


def ceilf(x):
    return -((-x.numerator) // x.denominator)


def up(x):
    return F(ceilf(x * GRID), GRID)


def pair(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def unpair(x):
    return F(int(x['numerator']), int(x['denominator']))


def afrac(x):
    assert x.is_exact() and x.is_finite()
    m, e = x.man_exp()
    m, e = int(m), int(e)
    return F(m * (1 << e)) if e >= 0 else F(m, 1 << (-e))


def fromdyad(x):
    m, e = int(x['mantissa']), int(x['exponent'])
    return F(m * (1 << e)) if e >= 0 else F(m, 1 << (-e))


def aball(x):
    return arb(x.numerator) / x.denominator


def ceilball(x):
    return int(x.upper().ceil().fmpz())


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while b := f.read(1024 * 1024):
            h.update(b)
    return h.hexdigest()


def primes(limit):
    return [q for q in range(2, limit + 1)
            if all(q % d for d in range(2, isqrt(q) + 1))]


def initialize():
    """Independent 384-bit roots, log enclosures, finite-exponent proof."""
    ctx.prec = 384
    active, records = [], []
    eps, sig = F(4, 45), F(1, 18)
    for q in primes(512):
        w = arb(q).root(45) ** (-4)
        r = arb(q).root(18) ** (-1)
        if (w * (1 + r)).upper() < 1:
            omitted = q
            omission = pair(afrac((w * (1 + r)).upper()))
            break
        E = 1
        while not (w ** E / (1 - r)).upper() < 1:
            E += 1
        wp = [arb(1)]
        rp, s0, s1, K, L = arb(1), arb(1), arb(0), arb(1), arb(0)
        for e in range(1, E):
            wp.append(wp[-1] * w)
            rp *= r
            s0 += rp
            s1 += e * rp
            K = max(K, (wp[-1] * s0).upper())
            L = max(L, (wp[-1] * s1).upper())
        logq = arb(q).log()
        active.append({'q': q, 'w': w, 'r': r, 'wp': wp, 'logq': logq})
        rec = {'q': q, 'first_discarded_exponent': E}
        for key, x in [('w', w), ('r', r), ('logq', logq)]:
            rec[key + '_lower'] = pair(afrac(x.lower()))
            rec[key + '_upper'] = pair(afrac(x.upper()))
        rec.update({'K_upper': pair(afrac(K)), 'L_upper': pair(afrac(L)),
                    'discarded_tail_upper': pair(afrac((w ** E / (1 - r)).upper()))})
        records.append(rec)
    assert omitted == 431 and len(active) == 82 and active[-1]['q'] == 421
    return active, {'epsilon': str(eps), 'sigma': str(sig),
        'local_precision': 384, 'active_primes': records,
        'first_omitted_prime': omitted, 'first_omitted_test_upper': omission,
        'pi_lower': pair(afrac(arb.pi().lower())),
        'pi_upper': pair(afrac(arb.pi().upper()))}


def read_complete(path, columns, panels, blocks, key, expected_raw):
    codec = struct.Struct('<' + 'Q' * columns)
    totals, bs, bi = [0] * columns, [0] * columns, 0
    h = hashlib.sha256()
    with gzip.open(path, 'rb') as stream:
        for node in range(panels):
            b = stream.read(codec.size)
            assert len(b) == codec.size
            h.update(b)
            values = codec.unpack(b)
            for j, value in enumerate(values):
                totals[j] += value
                bs[j] += value
            if node == blocks[bi]['last_node']:
                assert bs == blocks[bi][key]
                bs = [0] * columns
                bi += 1
        assert stream.read(1) == b'' and bi == len(blocks)
    assert h.hexdigest() == expected_raw
    return totals


def legacy_rank_bounds(parameters):
    rankpath = BASE / 'audit/interval_antichain/arb_178_full32768.json'
    suppath = BASE / 'audit/disjoint_support/support_169_32768.json'
    assert sha(rankpath) == EXPECTED_RANK
    assert sha(suppath) == EXPECTED_SUPPORT
    old = json.loads(rankpath.read_text())
    sup = json.loads(suppath.read_text())
    ri, si = old['integration'], sup['integration']
    assert si['source_sha256'] == EXPECTED_RANK
    assert old['epsilon'] == sup['epsilon'] == '4/45'
    assert old['sigma'] == sup['sigma'] == '1/18'
    q = [r['q'] for r in parameters['active_primes']]
    assert q == [r['q'] for r in old['parameters']['active_primes']]
    assert si['small_primes'] == [3, 5, 7, 11, 13, 17]
    assert si['large_tail_first_prime'] == 19
    rp, sp = rankpath.parent / ri['node_file'], suppath.parent / si['node_file']
    assert sha(rp) == ri['node_gzip_sha256']
    assert sha(sp) == si['gzip_sha256']
    N, scale, start, count = ri['panels'], 1 << ri['node_scale_bits'], 7, 64
    assert N == si['panels'] == 32768 and ri['node_scale_bits'] == si['node_scale_bits']
    rs = read_complete(rp, len(q) + 1, N, ri['blocks'], 'tail_sums', ri['node_raw_sha256'])
    ss = read_complete(sp, count, N, si['blocks'], 'subset_sums', si['raw_sha256'])
    assert rs == [int(t['node_sum_integer']) for t in ri['tails']]
    assert ss == [int(t['node_sum_integer']) for t in si['subset_integrals']]
    assert ri['node_raw_sha256'] == si['source_node_sha256']
    K = [unpair(r['K_upper']) for r in parameters['active_primes']]
    L = [unpair(r['L_upper']) for r in parameters['active_primes']]
    Delta = unpair(parameters['pi_upper']) / (2 * N)
    KK, FF, WW = [F(1)] * (len(q) + 1), [F(1)] * (len(q) + 1), [F(1)] * (len(q) + 1)
    for j in range(len(q) - 1, -1, -1):
        KK[j] = up(K[j] * KK[j + 1])
        FF[j] = up((1 + L[j] * Delta) * FF[j + 1])
        WW[j] = min(KK[j], up(F(rs[j], N * scale) * FF[j]))
    Ws, Ks = [], []
    for mask in range(count):
        cap, correction = KK[start], FF[start]
        for j in range(1, start):
            if mask & (1 << (j - 1)):
                cap = up(cap * K[j])
                correction = up(correction * (1 + L[j] * Delta))
        Ks.append(cap)
        Ws.append(min(cap, up(F(ss[mask], N * scale) * correction)))
    return {'rank_json': str(rankpath.relative_to(BASE)), 'support_json': str(suppath.relative_to(BASE)),
        'rank_json_sha256': EXPECTED_RANK, 'support_json_sha256': EXPECTED_SUPPORT,
        'rank_node_gzip_sha256': sha(rp), 'support_node_gzip_sha256': sha(sp),
        'rank_node_raw_sha256': ri['node_raw_sha256'], 'support_node_raw_sha256': si['raw_sha256'],
        'rank_values_checked': N * (len(q) + 1), 'support_values_checked': N * count,
        'rank_blocks_checked': len(ri['blocks']), 'support_blocks_checked': len(si['blocks']),
        'all_previous_nodes_and_block_sums_verified': True,
        'rank_K_upper': [pair(x) for x in KK], 'rank_W_upper': [pair(x) for x in WW],
        'subset_K_upper': [pair(x) for x in Ks], 'subset_W_upper': [pair(x) for x in Ws]}


def integrate(active, parameters, density, lambdas, path, blocksize):
    n, start, masks = len(active), 7, 64
    scale, columns = 1 << 24, n + 1 + masks
    N = max(lambdas) * density
    codec = struct.Struct('<' + 'Q' * columns)
    nodepath = path.with_suffix('.nodes.u64.gz')
    totals, block = [0] * columns, [0] * columns
    blocks, snapshots = [], []
    digest = hashlib.sha256()
    ctx.prec = 192
    started = time.monotonic()
    endpoints = {lam * density: lam for lam in lambdas}
    with gzip.open(nodepath, 'wb', compresslevel=6) as out:
        for node in range(N):
            t = arb(2 * node + 1) / (2 * density)
            gs = []
            for row in active:
                z = acb(0, t * row['logq']).exp()
                rz = row['r'] * z
                term, poly, g = acb(1), acb(1), arb(1)
                for e in range(1, len(row['wp'])):
                    term *= rz
                    poly += term
                    g = max(g, (row['wp'][e] * poly.abs_upper()).upper())
                gs.append(g)
            tails = [arb(1)] * (n + 1)
            values = [scale] * columns
            for j in range(n - 1, -1, -1):
                tails[j] = (tails[j + 1] * gs[j]).upper()
                values[j] = ceilball(tails[j] * scale)
            subsets = [tails[start]] * masks
            values[n + 1] = values[start]
            for mask in range(1, masks):
                low = mask & -mask
                bit = low.bit_length() - 1
                subsets[mask] = (subsets[mask ^ low] * gs[bit + 1]).upper()
                values[n + 1 + mask] = ceilball(subsets[mask] * scale)
            assert all(scale <= v < 2 ** 64 for v in values)
            assert values[n] == scale
            assert all(values[j] >= values[j + 1] for j in range(n))
            encoded = codec.pack(*values)
            out.write(encoded)
            digest.update(encoded)
            for j, v in enumerate(values):
                totals[j] += v
                block[j] += v
            if node + 1 in endpoints:
                snapshots.append({'lambda': endpoints[node + 1], 'panels': node + 1,
                                  'node_sums': totals.copy()})
            if (node + 1) % blocksize == 0 or node + 1 == N:
                blocks.append({'last_node': node, 'column_sums': block})
                block = [0] * columns
                progress = {'stage': 'Mellin all tails and supports', 'completed_nodes': node + 1,
                    'total_nodes': N, 'seconds_display_only': round(time.monotonic() - started, 3)}
                print(json.dumps(progress), flush=True)
                path.with_suffix('.checkpoint.json').write_text(json.dumps(progress) + '\n')
    # Re-read all node values, every block, every lambda prefix and the raw hash.
    rb, bs, bi, prefix_i = [0] * columns, [0] * columns, 0, 0
    readhash = hashlib.sha256()
    with gzip.open(nodepath, 'rb') as inp:
        for node in range(N):
            b = inp.read(codec.size)
            assert len(b) == codec.size
            readhash.update(b)
            values = codec.unpack(b)
            for j, v in enumerate(values):
                rb[j] += v
                bs[j] += v
            if node == blocks[bi]['last_node']:
                assert bs == blocks[bi]['column_sums']
                bs = [0] * columns
                bi += 1
            if node + 1 == snapshots[prefix_i]['panels']:
                assert rb == snapshots[prefix_i]['node_sums']
                prefix_i += 1
                if prefix_i == len(snapshots):
                    assert node + 1 == N
        assert inp.read(1) == b''
    assert rb == totals and readhash.digest() == digest.digest()
    assert bi == len(blocks) and prefix_i == len(snapshots)
    # Integral bounds have exact rational length and an outward rational correction.
    local = parameters['active_primes']
    K = [unpair(r['K_upper']) for r in local]
    D = [unpair(r['L_upper']) * unpair(r['logq_upper']) / (2 * density) for r in local]
    KK, FF = [F(1)] * (n + 1), [F(1)] * (n + 1)
    for j in range(n - 1, -1, -1):
        KK[j] = up(K[j] * KK[j + 1])
        FF[j] = up((1 + D[j]) * FF[j + 1])
    caps, corrections = KK.copy(), FF.copy()
    for mask in range(masks):
        cap, cor = KK[start], FF[start]
        for j in range(1, start):
            if mask & (1 << (j - 1)):
                cap = up(cap * K[j])
                cor = up(cor * (1 + D[j]))
        caps.append(cap)
        corrections.append(cor)
    for snap in snapshots:
        lam = snap['lambda']
        snap['integral_upper'] = [pair(min(up(lam * caps[j]),
            up(F(value, density * scale) * corrections[j])))
            for j, value in enumerate(snap['node_sums'])]
    return {'density': density, 'step': str(F(1, density)), 'panels': N,
        'max_lambda': max(lambdas), 'lambdas': lambdas, 'node_precision': 192,
        'columns': columns, 'tail_columns': n + 1, 'subset_columns': masks,
        'column_order': '82 prime tails in increasing first prime, empty tail, then masks 0..63',
        'small_primes': [3, 5, 7, 11, 13, 17], 'large_tail_first_prime': 19,
        'node_scale_bits': 24, 'node_file': nodepath.name, 'row_bytes': codec.size,
        'node_format': 'gzip little-endian uint64', 'raw_bytes': N * codec.size,
        'raw_sha256': digest.hexdigest(), 'gzip_sha256': sha(nodepath),
        'all_nodes_readback': True, 'values_checked': N * columns,
        'all_prefixes_readback': True, 'correction_upper': [pair(x) for x in corrections],
        'K_upper': [pair(x) for x in caps], 'prefix_integrals': snapshots,
        'blocks': blocks, 'seconds_display_only': time.monotonic() - started}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--density', type=int, default=16384)
    ap.add_argument('--lambdas', default='1,2,3,4,6,8,12,18')
    ap.add_argument('--block-size', type=int, default=4096)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    lambdas = sorted(set(map(int, args.lambdas.split(','))))
    assert args.density > 0 and lambdas and min(lambdas) > 0
    active, parameters = initialize()
    legacy = legacy_rank_bounds(parameters)
    print(json.dumps({'stage': 'setup', 'active_primes': len(active),
        'local_candidate_count': sum(len(r['wp']) - 1 for r in active),
        'old_rank_and_support_nodes_verified': True}), flush=True)
    integration = integrate(active, parameters, args.density, lambdas, path, args.block_size)
    result = {'status': 'INTEGRALS CERTIFIED; analytic Selberg theorem audited separately',
        'flint_version': flint.__version__, 'parameters': parameters,
        'legacy_rank_certificate': legacy, 'integration': integration,
        'script_sha256': sha(Path(__file__))}
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'panels': integration['panels'],
        'values_checked': integration['values_checked'], 'gzip_bytes': path.with_suffix('.nodes.u64.gz').stat().st_size,
        'seconds_display_only': integration['seconds_display_only']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
