"""Rebuild the full residual interval (1304, 10**66] and check all witnesses.

Uses only the Python standard library. Does not recompute large determinants.
The archived independent review supplies source binding and witness routines;
the new enumeration starts at 1304, not the inherited 10**56 endpoint.
"""
from pathlib import Path
import argparse, runpy, json, hashlib, math, time
from collections import Counter

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.output.exists():
        raise SystemExit('Choose a new output path; no overwrite.')
    root = a.project_root.resolve()
    reviewer = Path(__file__).resolve().parent / 'finite_review_portable.py'
    ns = runpy.run_path(str(reviewer), run_name='archive_functions')
    old = root / 'agents/finite_extension_round15'
    new = root / 'agents/finite_extension_round20'
    ns['sha'].__globals__.update(ROOT=root, OLD=old, NEW=new)
    start = time.monotonic()
    data = ns['load'](new / 'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json')
    ds = data['selected_d']
    assert len(ds) == len(set(ds)) == 57
    parameter_index = json.loads((reviewer.parent.parent / 'indices/general_rank_inputs.json').read_text(encoding='utf8'))
    assert ds == [entry['d'] for entry in parameter_index]
    rebuilt = ns['rebuild'](ds, 1304, 10**66)
    vals = rebuilt.pop('values')
    assert list(map(str, vals)) == data['candidate_list'] and len(vals) == 27985
    assert [int(row['n']) for row in data['records']] == vals
    assert all(n > 10**56 for n in vals)
    sequence_sha = hashlib.sha256(('\n'.join(map(str, vals))+'\n').encode()).hexdigest()
    assert sequence_sha == '5399b99b3aace9bf17fcf322d40dc8eaf3e2b8a2cb2cd7c9eb5a3c5afa4c9809'
    expected = {}
    composite_count = 0
    for row in data['records']:
        n = int(row['n'])
        if row['status'] == 'EXACT_COMPOSITE_WITNESS':
            ns['composite'](n, row['certificate'])
            composite_count += 1
        else:
            eligible = [d for d in (167,173) if n > 8*d and n % (2*d) not in (1,2*d-1)]
            assert eligible == row['eligible_d']
            expected[n] = eligible[0]
    assert composite_count == 26210
    assert Counter(expected.values()) == {167:1755,173:20}
    old_data = ns['load'](old / 'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json')
    reused = {int(n):expected[int(n)] for n in old_data['pointwise_candidates']}
    fresh = {n:d for n,d in expected.items() if n not in reused}
    sources = {d:ns['inspect_source'](d) for d in (167,173)}
    r1,m1,b1 = ns['pointwise'](old, ns['load'](old/'BATCH_PLAN.json'), old/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json', reused, sources)
    r2,m2,b2 = ns['pointwise'](new, ns['load'](new/'BATCH_PLAN.json'), new/'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json', fresh, sources)
    assert not m1 and not m2 and len(r1)+len(r2) == 1775
    assert all(r['status']=='IF_PRIME_THEN_CONNECTED' for r in [*r1.values(),*r2.values()])
    result = {'status':'PASS', 'lower_exclusive':1304, 'upper_inclusive':str(10**66),
        'crt':rebuilt, 'candidates':27985, 'composite_witnesses_recomputed':composite_count,
        'pointwise_source_bound_receipts':1775,
        'excluded_historical_executable_identities':list(ns['sha'].__globals__['SKIPPED_EXECUTABLE_IDENTITIES'].values()), 'batches':len(b1)+len(b2),
        'candidate_sequence_sha256':sequence_sha,
        'reviewer_sha256':hashlib.sha256(reviewer.read_bytes()).hexdigest(),
        'scope':'Fresh exhaustive CRT and modular compositeness checks, full source row parsing and saved unit determinant receipt checks; no fresh determinant reconstruction or general-parameter certification.',
        'seconds':time.monotonic()-start}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
