"""Recover complete integer coefficient rows from the 32 normalized matrices.

Recovery preserves the mathematical coefficients, not historical JSONL bytes.
The output is a candidate source for the separately supplied Phi recurrence;
the three low-moment identities alone do not certify source legality.
Uses only Python's standard library. Existing outputs are never overwritten.
"""
from pathlib import Path
from fractions import Fraction
import argparse, hashlib, json, math, sys, zipfile

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            r[i+j] += x*y
    return trim(r)


def divrem(p, q):
    p = [Fraction(x) for x in trim(p)]
    q = [Fraction(x) for x in trim(q)]
    a = [Fraction(0)] * max(1, len(p)-len(q)+1)
    while p != [0] and len(p) >= len(q):
        k = len(p)-len(q)
        c = p[-1]/q[-1]
        a[k] += c
        for i, x in enumerate(q):
            p[k+i] -= c*x
        p = trim(p)
    return trim(a), p


def monic(p):
    p = trim(p)
    return [x/p[-1] for x in map(Fraction, p)]


def filter_poly(d, t):
    """Rebuild the monic lcm filter from d,t, rather than trusting g_T."""
    primes, u = [], d
    for q in range(2, math.isqrt(d)+1):
        if u % q == 0:
            primes.append(q)
            while u % q == 0:
                u //= q
    if u > 1:
        primes.append(u)
    g = [Fraction(1)]
    for q in primes:
        def valuation(n):
            a = 0
            while n % q == 0:
                a += 1
                n //= q
            return a
        r = t // q**(valuation(t)-valuation(d)+1)
        n = r-1
        f = [0]*(n+1)
        for j in range(n//2+1):
            f[n-2*j] = (-1)**j * math.comb(n-j, j)
        if n % 2:
            f = [0]+f
        assert all(f[k] == 0 for k in range(1, len(f), 2))
        f = f[::2]
        a, b = g, f
        while b != [0]:
            _, rmd = divrem(a, b)
            a, b = b, rmd
        quotient, rmd = divrem(mul(g, f), monic(a))
        assert rmd == [0]
        g = monic(quotient)
    assert all(x.denominator == 1 for x in g)
    return [int(x) for x in g]


class Inputs:
    """Read and authenticate logical inputs from objects/ or the data ZIPs."""
    def __init__(self, package, input_tree=None):
        self.package, self.input_tree = package, input_tree
        self.rows = {x['path']:x for x in json.loads((package/'INPUTS.json').read_text(encoding='utf8'))}
        self.zips = []
        self.members = {}
        if input_tree is None:
            self.zips = [zipfile.ZipFile(p) for p in sorted((package/'data').glob('*.zip'))]
            for z in self.zips:
                for n in z.namelist():
                    if n in self.members:
                        raise ValueError('Duplicate object member: '+n)
                    self.members[n] = z

    def read(self, logical, expected=None):
        x = self.rows[logical]
        if expected is not None:
            assert x['sha256'] == expected, logical
        if self.input_tree is not None:
            b = (self.input_tree/logical).read_bytes()
        elif (self.package/'objects'/x['sha256']).is_file():
            b = (self.package/'objects'/x['sha256']).read_bytes()
        else:
            n = 'objects/'+x['sha256']
            b = self.members[n].read(n)
        assert len(b) == x['bytes'] and hashlib.sha256(b).hexdigest() == x['sha256'], logical
        return json.loads(b)


def restored_rows(entry, matrix):
    d, N, gate = entry['d'], entry['N'], entry['strict_p_gate']
    assert matrix['d'] == d and matrix['nd'] == N
    assert len(matrix['rows']) == len(matrix['row_metadata']) == entry['rows']
    parameters = {}
    for p in matrix['parameters']:
        t, g = p['n'], list(map(int, p['g_T']))
        assert t % d == 0 and d <= t <= N and t not in parameters
        assert g == filter_poly(d, t)
        assert p['m_count'] == t//2-(len(g)-1)
        parameters[t] = (g, p['m_count'])
    expected = {(t, m) for t, (_, count) in parameters.items() for m in range(count)}
    seen = set()
    for row, meta in zip(matrix['rows'], matrix['row_metadata']):
        t, m = int(meta['n']), int(meta['m'])
        assert (t, m) in expected and (t, m) not in seen
        seen.add((t, m))
        high = list(map(int, row))
        assert len(high) == N-2 and all(v == 0 for v in high[max(t-2, 0):])
        total = sum(math.comb(2*j, j)*v for j, v in enumerate(high, 3))
        assert total % 30 == 0
        c2 = -total//30
        normalized = [0, 12*c2, c2]+high
        content = int(meta['divided_small_content'])
        assert content > 0
        remainder = content
        for q in range(2, gate+1):
            while remainder % q == 0:
                remainder //= q
        assert remainder == 1, (d, t, m, 'content outside characteristic gate')
        full = [content*x for x in normalized[:t+1]]
        assert full[0] == 0 and full[1] == 12*full[2]
        assert sum(math.comb(2*j, j)*x for j, x in enumerate(full)) == 0
        yield {'d':d, 't':t, 'm':m, 'g_T':parameters[t][0],
               'divided_small_content':str(content), 'coefficients':[str(x) for x in full]}
    assert seen == expected


def source_height(t, g, m):
    # Absolute coefficient bound from f_t's explicit power-sum formula,
    # and l1(Phi(monomial)) <= 3^(total degree).
    s = 0
    for i in range(t+1):
        v = 2*t*math.comb(t+i, t-i)
        assert v % (t+i) == 0
        s += v//(t+i)*13**i
    return 9**(t+m)*s*sum(abs(x)*9**i for i, x in enumerate(g))


def write_jobs(path, blocks):
    """Produce exact-integer-comparison jobs for public phi_raw_modular.cpp.

    Every M is a power of two strictly greater than twice the proved height
    bound and the supplied coefficients. This is not random-prime sampling.
    """
    with path.open('x', encoding='utf8', newline='\n') as f:
        f.write(str(len(blocks))+'\n')
        for (d, t), rows in sorted(blocks.items()):
            g = rows[0]['g_T']
            assert all(r['g_T'] == g for r in rows)
            h = max(source_height(t, g, r['m']) for r in rows)
            M = 1 << (2*h).bit_length()
            assert M > 2*h
            f.write(f'{d} {t} {len(g)-1} {len(rows)}\n')
            f.write(' '.join(map(str, g))+'\n')
            for r in rows:
                c = list(map(int, r['coefficients']))
                assert max(map(abs, c)) <= source_height(t, g, r['m'])
                f.write(f"{r['m']} {M}\n")
                f.write(' '.join(str(x % M) for x in c)+'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--package', type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument('--index', type=Path)
    ap.add_argument('--input-tree', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--d', type=int, action='append', help='Select parameters; default is all 32.')
    ap.add_argument('--jobs', action='store_true', help='Also export complete Phi recurrence job streams.')
    a = ap.parse_args()
    index = a.index or a.package/'indices/RECOVERED_RANK_ROWS.json'
    entries = json.loads(index.read_text(encoding='utf8'))['entries']
    if a.d:
        assert set(a.d) <= {x['d'] for x in entries}
        entries = [x for x in entries if x['d'] in set(a.d)]
    if a.output.exists():
        raise SystemExit('Choose a new output directory; no overwrite.')
    a.output.mkdir(parents=True)
    reader = Inputs(a.package, a.input_tree)
    result = []
    for e in entries:
        cert = reader.read(e['certificate_path'], e['certificate_sha256'])
        assert cert['matrix_sha256'] == e['matrix_sha256']
        matrix = reader.read(e['matrix_path'], e['matrix_sha256'])
        assert matrix['raw_rows_sha256'] == e['historical_raw_sha256']
        assert matrix['raw_rows_file'] == Path(e['historical_raw_path']).name
        path = a.output/e['recovered_file']
        blocks, count, coefficients, h = {}, 0, 0, hashlib.sha256()
        with path.open('xb') as f:
            for row in restored_rows(e, matrix):
                b = (json.dumps(row, sort_keys=True, separators=(',', ':'))+'\n').encode()
                f.write(b);h.update(b);count += 1;coefficients += len(row['coefficients'])
                if a.jobs:
                    blocks.setdefault((row['d'], row['t']), []).append(row)
        if a.jobs:
            write_jobs(a.output/f"phi_d{e['d']}.in", blocks)
        item = {'d':e['d'], 'matrix_sha256':e['matrix_sha256'], 'rows':count,
                'full_integer_coefficients':coefficients, 'recovered_file':path.name,
                'recovered_sha256':h.hexdigest(), 'historical_raw_sha256':e['historical_raw_sha256'],
                'historical_bytes_reconstructed':False, 'source_recurrence_checked':False}
        result.append(item)
        print(json.dumps(item), flush=True)
    receipt = {'status':'PASS_COEFFICIENT_RECOVERY', 'entries':result,
               'scope':'Exact coefficient recovery and filter/parameter/content checks. Run the separate original Phi recurrence on exported jobs to certify source legality.'}
    (a.output/'RECOVERY_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')


if __name__ == '__main__':
    main()
