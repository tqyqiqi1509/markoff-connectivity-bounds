#!/usr/bin/env python3
"""Independent full-modulus signed-atom CRT enumeration and compositeness.

No old CRT producer is imported.  General criteria 139,149,151,157,163 remain
explicit gates; this program does not accept a finite interval.
"""
from pathlib import Path
import argparse, bisect, collections, datetime, hashlib, json, math, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LOWER = 10**56
BASES = (2, 3, 5, 7)
ADDITIONAL = (139, 149, 151, 157, 163)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temp.replace(path)


def factor_small(n):
    out, q = {}, 2
    while q*q <= n:
        while n % q == 0:
            out[q] = out.get(q, 0) + 1
            n //= q
        q += 1
    if n > 1:
        out[n] = 1
    return out


def components(ds):
    powers, edges = {}, []
    for d in ds:
        f = factor_small(2*d)
        for q, a in f.items():
            powers[q] = max(powers.get(q, 0), a)
        edges.append((d, [q for q, a in f.items() if q != 2 or a > 1]))
    graph = {q: set() for q in powers}
    for d, edge in edges:
        for q in edge:
            graph[q].update(set(edge) - {q})
    seen, parts = set(), []
    for q in sorted(graph):
        if q in seen:
            continue
        stack, part = [q], set()
        while stack:
            r = stack.pop()
            if r in part:
                continue
            part.add(r)
            stack.extend(graph[r] - part)
        seen |= part
        parts.append(sorted(part))
    blocks = [math.prod(q**powers[q] for q in part) for part in parts]
    M = math.prod(blocks)
    assert M == math.lcm(*(2*d for d in ds))
    assert all(math.gcd(x, y) == 1 for i, x in enumerate(blocks) for y in blocks[i+1:])
    assert powers[2] > 1, 'This campaign has a distinguishable 2-adic sign.'
    return M, powers, parts, blocks, edges


def signed_sums(atoms, M):
    out = [0]
    for atom in atoms:
        out = [(v+atom) % M for v in out] + [(v-atom) % M for v in out]
    assert len(out) == len(set(out)) == 2**len(atoms)
    return sorted(out)


def enumerate_interval(ds, lower, upper):
    M, powers, parts, blocks, edges = components(ds)
    assert 0 <= lower < upper < M
    atoms = [(M//b)*pow(M//b, -1, b) % M for b in blocks]
    assert sum(atoms) % M == 1
    for i, atom in enumerate(atoms):
        assert all(atom % b == int(i == j) for j, b in enumerate(blocks))
    # Contiguous halves of full-modulus signed atoms, independent from the
    # old reversed/interleaved half-modulus progression construction.
    split = len(atoms)//2
    left = signed_sums(atoms[:split], M)
    right = signed_sums(atoms[split:], M)
    values, row_counts, nxt = [], [], None
    for a in left:
        count = 0
        # 0 <= a+b < 2M, so exactly these two nonoverlapping lifts suffice.
        for turn in (0, 1):
            lo = bisect.bisect_right(right, turn*M+lower-a)
            hi = bisect.bisect_right(right, turn*M+upper-a)
            for b in right[lo:hi]:
                n = a+b-turn*M
                assert lower < n <= upper
                values.append(n)
            count += hi-lo
            pos = bisect.bisect_right(right, turn*M+upper-a)
            if pos < len(right):
                n = a+right[pos]-turn*M
                if upper < n < M:
                    nxt = n if nxt is None else min(nxt, n)
        row_counts.append(count)
    values.sort()
    assert len(values) == len(set(values)) == sum(row_counts)
    assert all(all(n % (2*d) in (1, 2*d-1) for d in ds) for n in values)
    sequence = ('\n'.join(map(str, values))+'\n').encode()
    return {
        'method': 'Full-modulus CRT idempotents, contiguous signed-atom halves, exact two-lift interval searches; no producer imports.',
        'selected_d': ds, 'lower_exclusive': str(lower), 'upper_inclusive': str(upper),
        'minimal_modulus': str(M), 'prime_exponents': powers,
        'sign_components': parts, 'block_moduli': list(map(str, blocks)),
        'coupling_edges': [{'d': d, 'active_sign_primes': edge} for d, edge in edges],
        'sign_component_count': len(parts), 'residue_class_count': len(left)*len(right),
        'half_counts': [len(left), len(right)], 'row_intersection_counts': row_counts,
        'candidate_count': len(values), 'candidate_list': list(map(str, values)),
        'candidate_sequence_sha256': hashlib.sha256(sequence).hexdigest(),
        'next_candidate_above_upper': str(nxt),
    }


def witness(n):
    assert n > 7 and n % 2
    d, s = n-1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for base in BASES:
        g = math.gcd(base, n)
        if 1 < g < n:
            return {'kind': 'proper_factor', 'factor': str(g), 'cofactor': str(n//g)}
        assert g == 1
        chain = [pow(base, d, n)]
        for _ in range(1, s):
            chain.append(chain[-1]*chain[-1] % n)
        if chain[0] != 1 and n-1 not in chain:
            return {'kind': 'strong_Miller_Rabin_failure', 'base': base,
                    'odd_exponent': str(d), 'two_adic_exponent': s,
                    'power_chain': list(map(str, chain))}
    return None


def classify(enum):
    counts, records, remaining = collections.Counter(), [], []
    assignments = collections.Counter()
    for text in enum['candidate_list']:
        n = int(text)
        w = witness(n)
        if w:
            r = {'n': text, 'status': 'EXACT_COMPOSITE_WITNESS', 'certificate': w}
        else:
            eligible = [d for d in (167, 173) if n > 8*d and n % (2*d) not in (1, 2*d-1)]
            r = {'n': text, 'status': 'POINTWISE_EVIDENCE_REQUIRED', 'primality_proved': False,
                 'preferred_d': eligible[0] if eligible else None, 'eligible_d': eligible}
            remaining.append(text)
            assignments[str(eligible[0]) if eligible else 'none'] += 1
        records.append(r)
        counts[r['status']] += 1
    return records, dict(counts), remaining, dict(assignments)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fallback', action='store_true', help='Prepare a separate conditional d163-free 10^63 backup; never replace the main target.')
    args = ap.parse_args()
    t = time.monotonic()
    accepted_path = ROOT/'control/FINITE_1E56_ROOT_ACCEPTANCE.json'
    accepted = json.loads(accepted_path.read_text())
    assert accepted['status'] == 'ACCEPTED_INTERNAL_PROOF_CHAIN'
    inherited_path = ROOT/'agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json'
    assert sha(inherited_path) == accepted['bindings'][str(inherited_path.relative_to(ROOT))]
    inherited = json.loads(inherited_path.read_text())
    base_ds = inherited['selected_d']
    assert len(base_ds) == 52 and 58 in base_ds and not set(base_ds) & set(ADDITIONAL)
    # Full independent calibration on the accepted 52-criterion list.
    calibration = enumerate_interval(base_ds, 1, LOWER)
    for key in ('candidate_list', 'candidate_sequence_sha256', 'minimal_modulus',
                'sign_components', 'residue_class_count', 'next_candidate_above_upper'):
        assert calibration[key] == inherited[key], ('CRT52 independent calibration', key)
    additional = [d for d in ADDITIONAL if not (args.fallback and d == 163)]
    ds = sorted(base_ds + additional)
    upper = 10**(63 if args.fallback else 66)
    result = enumerate_interval(ds, LOWER, upper)
    records, counts, remaining, assignments = classify(result)
    gates, pending = [], []
    for d in additional:
        p = ROOT/f'control/D{d}_ROOT_ACCEPTANCE.json'
        if p.exists():
            a = json.loads(p.read_text())
            assert a['status'] == 'ACCEPTED_GENERAL_PRIME_D_CRITERION' and a['d'] == d
            gates.append({'d': d, 'status': 'ROOT_ACCEPTANCE_EXISTS', 'path': str(p.relative_to(ROOT)), 'sha256': sha(p)})
        else:
            gates.append({'d': d, 'status': 'PENDING_ROOT_GENERAL_CRITERION'})
            pending.append(d)
    comparison = None
    if not args.fallback:
        p = ROOT/'planning/finite_extension/projected_crt57_1e66.json'
        old = json.loads(p.read_text())
        # Only compare after independent enumeration and compositeness work.
        assert result['candidate_list'] == old['candidate_list']
        assert result['candidate_sequence_sha256'] == old['candidate_sequence_sha256']
        assert result['next_candidate_above_upper'] == old['next_candidate_above_upper']
        assert remaining == old['untested_or_pseudoprime_candidates']
        comparison = {'planning_sha256': sha(p), 'candidate_list_agrees': True,
                      'next_boundary_agrees': True, 'unproved_survivors_agree': True}
    result.update({
        'status': 'CONDITIONAL_PENDING_CRITERIA_AND_POINTWISE_EVIDENCE',
        'interval_accepted': False, 'conjecture': 'OPEN',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'additional_general_gates': gates, 'pending_general_parameters_at_snapshot': pending,
        'full_crt52_calibration': {'count': calibration['candidate_count'], 'all_values_match': True,
                                 'sequence_sha256': calibration['candidate_sequence_sha256']},
        'classification_counts': counts, 'pointwise_count': len(remaining),
        'pointwise_source_assignment': assignments, 'pointwise_candidates': remaining,
        'records': records, 'planning_postcomparison': comparison,
        'inherited_finite_acceptance_sha256': sha(accepted_path),
        'inherited_crt52_sha256': sha(inherited_path), 'code_sha256': sha(__file__),
        'seconds': time.monotonic()-t,
        'scope': 'Exact arithmetic candidate enumeration and composite witnesses only. A four-base pass is never a primality assertion; every survivor needs a lawful determinant or further exact compositeness evidence. Uniform and pure-Vieta bridges remain explicit inherited dependencies.',
    })
    filename = 'FALLBACK_CRT56_1E63.json' if args.fallback else 'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json'
    save(HERE/filename, result)
    print(json.dumps({k: result[k] for k in ('status','candidate_count','classification_counts','pointwise_count','pointwise_source_assignment','pending_general_parameters_at_snapshot','seconds')}, indent=2))


if __name__ == '__main__':
    main()
