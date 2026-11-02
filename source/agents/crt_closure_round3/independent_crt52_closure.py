#!/usr/bin/env python3
"""Fresh CRT52 completeness and finite closure audit.

This file imports no producer.  It assembles CRT separately on two coprime
half-moduli, then intersects arithmetic progressions in the second modulus.
It consumes completed dual-determinant receipts, but performs no determinant.
The inherited uniform rank theorems and the old accepted finite interval are
explicit mathematical dependencies, not independently reproved by this code.
"""
from pathlib import Path
import argparse, bisect, collections, datetime, hashlib, json, math, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / 'prior_restored'
ASSESS = PRIOR / 'latest/finite51_interval_assessment'
FINITE = ROOT / 'agents/attachment_audit/finite52'
UPPER = 10**56
OLD_UPPER = 15 * 10**54


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda: stream.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(name, value):
    p = HERE / name
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(p)


def factors(n):
    result = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = 1
    return result


def sign_partition(ds):
    """Merge overlapping hyperedges as sets; do not use a DSU or graph walk."""
    modulus = math.lcm(*(2*d for d in ds))
    prime_exponents = factors(modulus)
    parts = [{p} for p in sorted(prime_exponents)]
    edges = []
    for d in ds:
        # Signs +1 and -1 agree modulo 2, but not modulo 4 or any odd prime.
        support = {p for p, e in factors(2*d).items() if p != 2 or e >= 2}
        edges.append({'d': d, 'active_sign_primes': sorted(support)})
        overlap = [part for part in parts if part & support]
        if overlap:
            union = set().union(*overlap)
            parts = [part for part in parts if not part & support] + [union]
    parts = sorted((sorted(part) for part in parts), key=lambda x: x[0])
    blocks = [math.prod(p**prime_exponents[p] for p in part) for part in parts]
    assert math.prod(blocks) == modulus
    assert all(math.gcd(a, b) == 1 for i, a in enumerate(blocks) for b in blocks[i+1:])
    # An independent exhaustive check of the partition equivalence is cheap:
    # each original sign edge is contained in a final part, and each merger
    # was made only because an actual criterion connected those parts.
    assert all(any(set(edge['active_sign_primes']) <= set(part) for part in parts)
               for edge in edges if edge['active_sign_primes'])
    return modulus, prime_exponents, parts, blocks, edges


def half_crt(blocks):
    """Incremental CRT in the small half-modulus, not full-modulus atoms."""
    modulus, residues = 1, [0]
    for q in blocks:
        inv = pow(modulus, -1, q)
        signs = sorted({1 % q, -1 % q})
        residues = [u + modulus * (((s-u)*inv) % q)
                    for u in residues for s in signs]
        modulus *= q
        assert all(0 <= u < modulus for u in residues)
        assert len(residues) == len(set(residues))
    return modulus, sorted(residues)


def enumerate_progressions(ds, upper):
    modulus, prime_exponents, parts, blocks, edges = sign_partition(ds)
    assert 2 <= upper < modulus
    # Reverse and alternate block order, distinct from both inherited splits.
    ordered = blocks[::-1]
    A, left = half_crt(ordered[::2])
    B, right = half_crt(ordered[1::2])
    assert A*B == modulus and math.gcd(A, B) == 1
    inverse = pow(A, -1, B)
    wlist = sorted({v*inverse % B for v in right})
    assert len(wlist) == len(right)
    values = []
    next_value = None
    row_counts = []
    for u in left:
        # n=u+A*t, 0<=t<B, and t=(w-u/A) modulo B.
        first_t = max(0, -((u-2)//A))
        last_t = min(B-1, (upper-u)//A)
        offset = u*inverse % B
        count = 0
        if first_t <= last_t:
            low, high = offset+first_t, offset+last_t
            for shift in (0, B):
                lo = bisect.bisect_left(wlist, low-shift)
                hi = bisect.bisect_right(wlist, high-shift)
                for w in wlist[lo:hi]:
                    t = w+shift-offset
                    n = u+A*t
                    assert 2 <= n <= upper and 0 <= t < B
                    values.append(n)
                count += hi-lo
        # First point of this arithmetic progression after the interval.
        t0 = max(0, (upper-u)//A + 1)
        if t0 < B:
            threshold = offset+t0
            for shift in (0, B):
                k = bisect.bisect_left(wlist, threshold-shift)
                if k < len(wlist):
                    t = wlist[k]+shift-offset
                    n = u+A*t
                    if t0 <= t < B and upper < n < modulus:
                        next_value = n if next_value is None else min(next_value, n)
        row_counts.append(count)
    values.sort()
    assert len(values) == len(set(values)) == sum(row_counts)
    assert all(all(n % (2*d) in (1, 2*d-1) for d in ds) for n in values)
    seqsha = hashlib.sha256(('\n'.join(map(str, values))+'\n').encode()).hexdigest()
    return {
        'selected_d': ds, 'minimal_modulus': str(modulus),
        'prime_exponents': prime_exponents, 'sign_components': parts,
        'block_moduli': list(map(str, blocks)), 'coupling_edges': edges,
        'sign_component_count': len(parts),
        'residue_class_count': len(left)*len(right),
        'half_moduli': [str(A), str(B)], 'half_residue_counts': [len(left), len(right)],
        'upper_inclusive': str(upper), 'candidate_count': len(values),
        'candidate_sequence_sha256': seqsha, 'candidate_list': list(map(str, values)),
        'next_candidate_above_upper': str(next_value), 'row_intersection_counts': row_counts,
        'method': 'Incremental CRT on coprime half-moduli; exact progression intervals n=u+A*t; cyclic right-residue search.',
    }


def composite_witness(n, witness):
    """Check failure witnesses by iterative squaring, never infer primality."""
    if witness['kind'] == 'proper_factor':
        f = int(witness['factor'])
        assert 1 < f < n and f*int(witness['cofactor']) == n
        return
    assert witness['kind'] == 'strong_Miller_Rabin_failure'
    a, d, s = int(witness['base']), int(witness['odd_exponent']), witness['two_adic_exponent']
    assert 1 < a < n and math.gcd(a, n) == 1
    assert d > 0 and d % 2 and s > 0 and d*2**s == n-1
    chain = [pow(a, d, n)]
    for _ in range(1, s):
        chain.append(chain[-1]**2 % n)
    assert chain == list(map(int, witness['power_chain']))
    assert chain[0] != 1 and n-1 not in chain


def collect_completed(plan, expected):
    """Read only atomic COMPLETE files and bind both original numeric receipts."""
    records, sources, bindings = {}, {}, []
    # Each input source is read/hash checked once, despite repeated symlinks.
    source_hashes = {
        139: 'db7719d99be914a193814d98fdeacc3589d76b33dfe65cb2fdcd6d78951d219d',
        149: '1fb4aa3ff2ca819695d8625814baee0a1d9bb8321424f9ba0a0b9f9cb85ee6e3',
    }
    for d in source_hashes:
        p = (FINITE/plan['pilot_directory']/f'd{d}_prefix_rows.jsonl' if d == 139
             else FINITE/'sources'/f'd{d}_prefix_rows.jsonl')
        assert sha(p) == source_hashes[d]
        sources[d] = load(PRIOR/'latest/finite'/f'd{d}_source_audit.json')
        assert sources[d]['source_sha256'] == source_hashes[d]
    groups = [(plan['pilot_directory'], 139, None)]
    groups += [(r['directory'], r['d'], r) for r in plan['batches']]
    for directory_name, d, batch in groups:
        directory = FINITE/directory_name
        result_path = directory/f'd{d}_complete_verification.json'
        if not result_path.exists():
            continue
        result = load(result_path)
        assert result['status'] == 'COMPLETE' and result['d'] == d and result['N'] == 4*d
        if batch is None:
            assert sha(result_path) == plan['pilot_complete_sha256']
            manifest_path = directory/'pilot_manifest.json'
        else:
            assert sha(directory/'batch_manifest.json') == batch['manifest_sha256']
            manifest_path = directory/'batch_manifest.json'
        manifest = load(manifest_path)
        assert len(manifest['groups']) == 1
        group = manifest['groups'][0]
        assert group['d'] == d and group['N'] == 4*d
        listed = list(map(int, group['candidates']))
        if batch is not None:
            assert listed == list(map(int, batch['candidates']))
        for name, digest in plan['common_code_sha256'].items():
            assert sha(directory/name) == digest
        assert result['gmp_code_sha256'] == plan['common_code_sha256']['independent_modular_det.cpp']
        assert result['python_code_sha256'] == plan['common_code_sha256']['run_manifest_group_verification.py']
        assert result['source_sha256'] == source_hashes[d]
        audit = load(directory/f'd{d}_source_audit.json')
        assert audit['source_sha256'] == source_hashes[d]
        assert audit['all_three_moment_relations'] and audit['all_rows_legal_and_complete']
        assert audit['source_row_count'] == 4*d-2 and audit['row_order'] == sources[d]['row_order']
        assert audit['candidate_manifest_sha256'] == sha(manifest_path)
        assert audit['generator_sha256'] == plan['common_code_sha256']['phi_stream_selected.cpp']
        cross = [json.loads(line) for line in (directory/f'd{d}_cross_verified_certificates.jsonl').read_text().splitlines()]
        gmp = [json.loads(line) for line in (directory/f'd{d}_gmp_determinants.jsonl').read_text().splitlines()]
        assert cross == result['records']
        assert [int(x['candidate']) for x in cross] == [int(x['candidate']) for x in gmp] == listed
        assert len(cross) == result['candidate_count'] == group['candidate_count']
        counts = collections.Counter()
        for record, g in zip(cross, gmp):
            assert all(record[k] == value for k, value in g.items())
            n = int(record['candidate'])
            assert n in expected and expected[n] == d and n not in records
            assert record['d'] == d and record['N'] == 4*d and record['source_sha256'] == source_hashes[d]
            assert n > 8*d and n % (2*d) not in (1, 2*d-1)
            assert record['characteristic_gate'] and record['congruence_gate']
            status = record['accepted_status']
            if status == 'IF_PRIME_THEN_CONNECTED':
                v = int(record['determinant_mod_n'])
                assert record['status'] == 'nonzero_determinant' and 0 < v < n and math.gcd(v, n) == 1
                assert v == int(record['flint_determinant_mod_n'])
                assert record['unit_pivots'] == 4*d-2 and record['independent_determinants_agree']
            elif status == 'EXACT_COMPOSITE_FACTOR':
                f = int(record['factor'])
                assert record['status'] == 'composite_factor' and 1 < f < n and n % f == 0
            else:
                assert status == 'UNRESOLVED_ZERO_SELECTED_MINOR'
                assert int(record['determinant_mod_n']) == int(record['flint_determinant_mod_n']) == 0
            counts[status] += 1
            records[n] = dict(record, complete_receipt_sha256=sha(result_path))
        assert result['conditional_connected'] == counts['IF_PRIME_THEN_CONNECTED']
        assert result['composite'] == counts['EXACT_COMPOSITE_FACTOR']
        assert result['unresolved'] == counts['UNRESOLVED_ZERO_SELECTED_MINOR']
        bindings.append({'path': str(result_path.relative_to(ROOT)), 'sha256': sha(result_path), 'count': len(cross)})
    return records, bindings


def main():
    start = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    old_enum_path = ASSESS/'crt51_enumeration_1e56.json'
    old_assess_path = ASSESS/'crt51_workload_assessment_1e56.json'
    old_closed_path = PRIOR/'latest/finite51_15e54/complete_interval_closure_15e54.json'
    baseline_path = PRIOR/'latest/root/FINAL_ACCEPTANCE_STATUS.json'
    baseline = load(baseline_path)
    binding = {r['path']: r for r in baseline['bindings']}
    for p in [old_enum_path, old_assess_path, old_closed_path]:
        b = binding[str(p.relative_to(PRIOR))]
        assert p.stat().st_size == b['size'] and sha(p) == b['sha256']
    old_enum, old_assess, old_closed = map(load, [old_enum_path, old_assess_path, old_closed_path])
    assert old_closed['interval_accepted'] and old_closed['remaining_count'] == 0
    assert int(old_closed['upper_exclusive']) == OLD_UPPER
    assert int(baseline['accepted_finite_range']['integer_upper_exclusive']) == OLD_UPPER
    ds51 = old_enum['selected_d']
    assert len(ds51) == 51 and 58 not in ds51 and not set(ds51) & {139, 149, 151}
    old50_path = ASSESS/'crt_50_interval_1e55.json'
    assert sha(old50_path) == old_enum['input_bindings']['crt50']
    old50 = load(old50_path)
    assert old50['status'] == 'ALL_50_RANK_INPUTS_ACCEPTED_AND_CRT_REBUILT'
    assert old50['largest_new_strict_characteristic_bound'] == 1096 < OLD_UPPER
    assert sorted(old50['selected_d']+[46]) == ds51
    d46_path = PRIOR/'latest/finite46_recomputed/composite_d46_independent_verification.json'
    assert sha(d46_path) == old_enum['input_bindings']['d46_independent']
    d46 = load(d46_path)
    assert d46['strict_p_bound'] == 552 < OLD_UPPER and d46['excluded_residues_modulus'] == 92
    rebuilt51 = enumerate_progressions(ds51, UPPER)
    for key in ['minimal_modulus', 'candidate_list', 'candidate_sequence_sha256', 'candidate_count', 'sign_components', 'next_candidate_above_upper']:
        assert rebuilt51[key] == old_enum[key], key
    d58_path = ROOT/'control/D58_ROOT_ACCEPTANCE.json'
    d58 = load(d58_path)
    assert d58['status'] == 'ACCEPTED_GENERAL_D58_CRITERION' and d58['d'] == 58 and d58['N'] == 348
    for path, digest in d58['input_bindings'].items():
        assert sha(ROOT/path) == digest
    current = enumerate_progressions(sorted(ds51+[58]), UPPER)
    filtered = [n for n in rebuilt51['candidate_list'] if int(n) % 116 in (1, 115)]
    assert current['candidate_list'] == filtered
    assert current['candidate_count'] == 6988 and current['sign_component_count'] == 24
    assert current['minimal_modulus'] == old_enum['minimal_modulus']
    assert current['residue_class_count'] == 2**24
    assert current['sign_components'][0] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    current.update({'status': 'CRT52_COMPLETE_INDEPENDENT_PROGRESSION_ENUMERATION',
                    'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    'old51_fresh_full_list_agrees': True, 'old51_count': rebuilt51['candidate_count'],
                    'uniform139_149_151_used': False,
                    'input_bindings': {str(p.relative_to(ROOT)): sha(p) for p in [old_enum_path, d58_path]},
                    'source_sha256': sha(__file__)})
    save('CRT52_INDEPENDENT_ENUMERATION.json', current)
    assessed = {int(r['n']): r for r in old_assess['records']}
    assert sorted(assessed) == list(map(int, rebuilt51['candidate_list']))
    old_by = {int(r['n']): r for r in old_closed['records']}
    assert sorted(old_by) == [int(n) for n in rebuilt51['candidate_list'] if int(n) < OLD_UPPER]
    expected = {}
    for s in current['candidate_list']:
        n = int(s)
        if n >= OLD_UPPER and assessed[n]['status'].startswith('PENDING'):
            legal = [d for d in (139, 149) if n > 8*d and n % (2*d) not in (1, 2*d-1)]
            assert legal
            expected[n] = legal[0]
    assert len(expected) == 401 and collections.Counter(expected.values()) == {139: 394, 149: 7}
    plan_path = FINITE/'BATCH_PLAN.json'
    plan = load(plan_path)
    actual, completed_bindings = collect_completed(plan, expected)
    counts, records, missing, unresolved = collections.Counter(), [], [], []
    new_witnesses = 0
    old_closed_sha = sha(old_closed_path)
    for s in current['candidate_list']:
        n = int(s)
        if n < OLD_UPPER:
            assert n in old_by
            record = {'n': s, 'status': 'INHERITED_ACCEPTED_15E54_RANGE',
                      'old_record_status': old_by[n]['status'], 'old_closure_sha256': old_closed_sha}
        elif assessed[n]['status'] == 'EXACT_COMPOSITE_WITH_NEW_WITNESS':
            witness = assessed[n]['certificate']
            composite_witness(n, witness)
            record = {'n': s, 'status': 'EXACT_COMPOSITE_RECHECKED', 'certificate': witness}
            new_witnesses += 1
        else:
            assert n in expected
            if n not in actual:
                record = {'n': s, 'status': 'MISSING_COMPLETED_POINTWISE_RECEIPT', 'd': expected[n]}
                missing.append(s)
            else:
                a = actual[n]
                status = a['accepted_status']
                record = {'n': s, 'status': status, 'd': a['d'], 'N': a['N'],
                          'complete_receipt_sha256': a['complete_receipt_sha256']}
                if status == 'IF_PRIME_THEN_CONNECTED':
                    record.update({'determinant_mod_n': a['determinant_mod_n'],
                                   'source_sha256': a['source_sha256'], 'primality_proved': False})
                elif status == 'EXACT_COMPOSITE_FACTOR':
                    record['factor'] = a['factor']
                else:
                    unresolved.append(s)
        counts[record['status']] += 1
        records.append(record)
    assert len(records) == len({r['n'] for r in records}) == 6988
    complete = not missing and not unresolved
    output = {
        'status': 'ALL_CRT52_NECESSARY_INTEGERS_THROUGH_1E56_CLOSED' if complete else 'PARTIAL_CRT52_CLOSURE_NO_INTERVAL_PROMOTION',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'interval_accepted': complete, 'upper_inclusive': str(UPPER),
        'continuous_prime_conclusion': 'Every prime 5<=p<=10^56 has connected pure-three-Vieta Markoff graph.' if complete else None,
        'inherited_finite_upper_exclusive': str(OLD_UPPER),
        'criterion_count': 52, 'new_uniform_parameter': 58, 'uniform139_149_151_used': False,
        'candidate_count': len(records), 'counts': dict(counts),
        'new_exact_composite_witnesses_recomputed': new_witnesses,
        'expected_pointwise': 401, 'completed_pointwise': len(actual),
        'missing_count': len(missing), 'unresolved_count': len(unresolved),
        'missing': missing, 'unresolved': unresolved,
        'source_sha256': sha(__file__),
        'input_bindings': {str(p.relative_to(ROOT)): sha(p) for p in [old_enum_path, old_assess_path, old_closed_path, baseline_path, old50_path, d46_path, d58_path, plan_path]},
        'completed_receipts': completed_bindings,
        'enumeration_sha256': sha(HERE/'CRT52_INDEPENDENT_ENUMERATION.json'),
        'scope': 'Fresh CRT reconstruction and exact composite witness rechecks. The inherited small-prime base and uniform rank-to-pure-Vieta theorems remain explicit dependencies. Completed determinant receipts are bound and their consistency/gates checked, not numerically recomputed here. No probable-prime pass is accepted as primality.',
        'seconds': time.monotonic()-start, 'records': records,
    }
    save('CRT52_FINITE_CLOSURE_REVIEW.json', output)
    print(json.dumps({k: v for k, v in output.items() if k not in ['records', 'missing', 'unresolved', 'completed_receipts', 'input_bindings']}, indent=2), flush=True)
    if args.require_complete:
        assert complete, 'Missing or unresolved completed pointwise evidence; no interval promotion.'


if __name__ == '__main__':
    main()
