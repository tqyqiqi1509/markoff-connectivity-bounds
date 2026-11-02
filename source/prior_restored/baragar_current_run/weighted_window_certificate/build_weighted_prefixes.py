#!/usr/bin/env python3
"""Exact prefixes for weighted exponential-window min(A,B/t) certification.

Reads the existing complete Mellin node stream. The weight in cell i>=1 is
bounded at its left endpoint i/density; i=0 is bounded separately by L.
Sums of V_i/i are bounded in both directions on a 2^-128 grid, so differences
use an upper final prefix minus a lower initial prefix. No upper-minus-upper
subtraction and no re-evaluation/interpolation of analytic nodes is used.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import struct
import sys
import time

sys.set_int_max_str_digits(0)
BITS = 128
S = 1 << BITS


def rat(x):
    return F(int(x['numerator']), int(x['denominator']))


def pair(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator)}


def ceilf(x):
    return -((-x.numerator) // x.denominator)


def up(x):
    return F(ceilf(x * S), S)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        while b := stream.read(1024 * 1024):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    path, out = Path(args.input), Path(args.output)
    d = json.loads(path.read_text())
    ps, info = d['parameters']['active_primes'], d['integration']
    density, N, columns = info['density'], info['panels'], info['columns']
    thresholds, universal, needed = [], [], {r['panels'] for r in info['prefix_integrals']}
    K = [rat(r['K_upper']) for r in ps]
    for j, row in enumerate(ps):
        cap = F(1)
        for k in range(j + 1, len(ps)):
            cap = up(cap * K[k])
        # An intentionally loose common cap covers small rounding-order changes
        # in both independent recurrences; every actual tail is checked below it.
        cap = 2 * (cap + 1)
        w = rat(row['w_upper'])
        A = 1
        while w ** A * (cap + 1) > 1:
            A += 1
        universal.append({'q': row['q'], 'first_omitted_exponent': A,
                          'Euler_tail_cap': pair(cap)})
        sigma = F(d['parameters']['sigma'])
        rlo, rhi = rat(row['r_lower']), rat(row['r_upper'])
        rholo, rhohi = F(1), F(1)
        for a in range(1, A):
            val = rholo*rlo
            rholo = F((val.numerator*S)//val.denominator,S)
            rhohi = up(rhohi*rhi)
            L = (1-rholo)/sigma
            B = 1+rhohi
            cut = B*density/L
            split = cut.numerator//cut.denominator+1
            assert split >= 1 and 0 < rholo <= rhohi < 1
            thresholds.append({'q': row['q'], 'a': a, 'split_panels': split,
                               'L_upper': pair(L), 'B_upper': pair(B),
                               'rho_lower':pair(rholo), 'rho_upper':pair(rhohi)})
            needed.add(min(split, N))
    needed = sorted(needed)
    nodepath = path.parent / info['node_file']
    assert sha(nodepath) == info['gzip_sha256']
    codec = struct.Struct('<' + 'Q' * columns)
    assert codec.size == info['row_bytes']
    sums, lows, highs, block = [0] * columns, [0] * columns, [0] * columns, [0] * columns
    prefixes, digest, bi, ni = [], hashlib.sha256(), 0, 0
    started = time.monotonic()
    with gzip.open(nodepath, 'rb') as stream:
        for i in range(N):
            data = stream.read(codec.size)
            assert len(data) == codec.size
            digest.update(data)
            values = codec.unpack(data)
            for j, v in enumerate(values):
                sums[j] += v
                block[j] += v
                if i:
                    lo, rem = divmod(v << BITS, i)
                    lows[j] += lo
                    highs[j] += lo + bool(rem)
            if i + 1 == needed[ni]:
                prefixes.append({'panels': i + 1, 'node_sums': list(map(str, sums)),
                    'harmonic_lower_scaled': list(map(str, lows)),
                    'harmonic_upper_scaled': list(map(str, highs))})
                ni += 1
                if ni == len(needed):
                    assert i + 1 == N
            if i == info['blocks'][bi]['last_node']:
                assert block == info['blocks'][bi]['column_sums']
                block = [0] * columns
                bi += 1
            if (i + 1) % 32768 == 0 or i + 1 == N:
                print(json.dumps({'stage': 'exact harmonic prefixes', 'completed_nodes': i + 1,
                    'total_nodes': N, 'seconds_display_only': round(time.monotonic() - started, 3)}), flush=True)
        assert stream.read(1) == b''
    assert digest.hexdigest() == info['raw_sha256'] and bi == len(info['blocks'])
    assert ni == len(needed)
    lookup = {r['panels']: r for r in prefixes}
    for pref in info['prefix_integrals']:
        assert list(map(int, lookup[pref['panels']]['node_sums'])) == pref['node_sums']
    prefixpath = out.with_suffix('.prefixes.json.gz')
    raw = json.dumps(prefixes, separators=(',', ':')).encode()
    with gzip.open(prefixpath, 'wb', compresslevel=6) as stream:
        stream.write(raw)
    with gzip.open(prefixpath, 'rb') as stream:
        rb = stream.read()
    assert rb == raw and json.loads(rb) == prefixes
    meta = {'status': 'ALL WEIGHT PREFIXES CERTIFIED',
        'input_file': str(path), 'input_sha256': sha(path), 'input_node_raw_sha256': info['raw_sha256'],
        'input_node_gzip_sha256': info['gzip_sha256'], 'source_values_read': N * columns,
        'all_source_blocks_verified': True, 'all_source_lambda_sums_verified': True,
        'harmonic_scale_bits': BITS, 'density': density, 'columns': columns, 'panels': N,
        'prefix_count': len(prefixes), 'thresholds': thresholds, 'universal_exponent_bounds': universal,
        'prefix_file': prefixpath.name, 'prefix_gzip_sha256': sha(prefixpath),
        'prefix_raw_sha256': hashlib.sha256(raw).hexdigest(), 'prefix_raw_bytes': len(raw),
        'all_prefix_records_readback': True, 'script_sha256': sha(Path(__file__)),
        'seconds_display_only': time.monotonic() - started}
    out.write_text(json.dumps(meta, indent=2) + '\n')
    print(json.dumps({k: v for k, v in meta.items() if k not in ['thresholds', 'universal_exponent_bounds']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
