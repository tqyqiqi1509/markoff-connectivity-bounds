#!/usr/bin/env python3
"""Independent arithmetic/source/receipt closure, without new determinants.

Rebuilds CRT using the older half-modulus progression implementation, distinct
from this branch's full-modulus signed-atom producer. Every composite power is
recomputed directly rather than by the producer's iterative squaring chain.
All source coefficients and moments are parsed/rechecked once per ready d.
Root remains responsible for final mathematical acceptance.
"""
from pathlib import Path
import argparse, collections, datetime, hashlib, importlib.util, json, math, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json'
PLAN = HERE/'BATCH_PLAN.json'
ADDITIONAL = (139,149,151,157,163)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path,obj):
    p=path.with_suffix(path.suffix+'.tmp')
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
    p.replace(path)


def check_composite(n,w):
    if w['kind']=='proper_factor':
        f=int(w['factor'])
        assert 1<f<n and f*int(w['cofactor'])==n
        return
    assert w['kind']=='strong_Miller_Rabin_failure'
    a,d,s=int(w['base']),int(w['odd_exponent']),w['two_adic_exponent']
    assert 1<a<n and math.gcd(a,n)==1 and d>0 and d%2 and s>0
    assert n-1==d*2**s
    # Direct independent exponents, not a recurrence through prior residues.
    actual=[pow(a,d*2**j,n) for j in range(s)]
    assert actual==list(map(int,w['power_chain']))
    assert actual[0]!=1 and n-1 not in actual


def source_check(d):
    ready_path=HERE/f'd{d}_SOURCE_READY.json'
    if not ready_path.exists():
        return None
    ready=load(ready_path)
    assert ready['d']==d and ready['N']==4*d
    assert ready['actual_generator']=='phi_stream_partitioned.cpp'
    for name,digest in ready['bindings'].items():
        assert sha(ROOT/name)==digest
    directory=ROOT/'new_pointwise'/f'd{d}'
    assert sha(directory/'phi_stream_partitioned.cpp')==ready['generator_sha256']
    assert sha(directory/'phi_stream_selected.cpp')==ready['legacy_selected_generator_equivalent_prefix_sha256']
    source=ROOT/ready['source_path']
    assert source.stat().st_size==ready['source_bytes'] and sha(source)==ready['source_sha256']
    N=4*d
    expected=[(b*d,m) for b in range(1,5) for m in range(b*(d-1)//2 if b<4 else d+1)]
    assert len(expected)==N-2
    expected_set=set(expected)
    binom=[math.comb(2*j,j) for j in range(N+1)]
    metadata,coefficient_count={},0
    with source.open() as f:
        for lineno,line in enumerate(f):
            assert line.endswith('\n')
            row=json.loads(line)
            assert row['d']==d and row['nd']==N
            key=(row['n'],row['m'])
            assert key in expected_set and key not in metadata
            coeff=row['coefficients']
            assert len(coeff)==row['n']+1
            assert all(isinstance(v,str) and str(int(v))==v for v in coeff)
            v=list(map(int,coeff))
            assert v[0]==0 and v[1]==12*v[2]
            assert sum(v[j]*binom[j] for j in range(len(v)))==0
            assert row['max_half_total_degree']==2*row['n']+(row['n']//d)//2+row['m']
            metadata[key]={'source_line_zero_based':lineno,'n':key[0],'m':key[1],
                'full_coefficients_sha256':hashlib.sha256((','.join(coeff)+'\n').encode()).hexdigest()}
            coefficient_count+=len(v)
    assert set(metadata)==expected_set
    ordered=[dict(metadata[key],matrix_row=i) for i,key in enumerate(expected)]
    result={'status':'ALL_STORED_INTEGER_COEFFICIENTS_MOMENTS_AND_ROW_KEYS_RECHECKED',
            'd':d,'N':N,'source_sha256':ready['source_sha256'],'source_bytes':ready['source_bytes'],
            'row_count':len(ordered),'coefficient_count':coefficient_count,'row_order':ordered,
            'source_ready_sha256':sha(ready_path),'scope':'Fresh complete numeric parsing and moment/row-order checks; not a second polynomial generation or determinant run.'}
    save(HERE/f'd{d}_INDEPENDENT_SOURCE_AUDIT.json',result)
    return result


def read_lines(path):
    result=[]
    with path.open() as f:
        for line in f:
            assert line.endswith('\n')
            result.append(json.loads(line))
    return result


def completed_pointwise(plan,expected,sources):
    actual,receipts={},[]
    for batch in plan['batches']:
        d=batch['d'];directory=HERE/batch['directory']
        path=directory/f'd{d}_complete_verification.json'
        if not path.exists():
            continue
        assert sources[d] is not None
        source=sources[d]
        manifest_path=directory/'batch_manifest.json'
        assert sha(manifest_path)==batch['manifest_sha256']
        manifest=load(manifest_path)
        assert len(manifest['groups'])==1
        group=manifest['groups'][0]
        assert group['d']==d and group['N']==4*d
        assert group['candidates']==batch['candidates']
        assert group['candidate_count']==len(group['candidates'])==batch['candidate_count']
        assert manifest['classification_sha256']==sha(DATA)
        assert (directory/f'd{d}_candidates.txt').read_text().splitlines()==batch['candidates']
        for name,digest in plan['common_code_sha256'].items():
            assert sha(directory/name)==digest
        # The actual Python file location matters because its BASE uses resolve.
        wrapper=directory/'run_manifest_group_verification.py'
        assert wrapper.resolve().parent==directory.resolve()
        result=load(path)
        assert result['status']=='COMPLETE' and result['d']==d and result['N']==4*d
        assert result['source_sha256']==source['source_sha256']
        assert result['gmp_code_sha256']==plan['common_code_sha256']['independent_modular_det.cpp']
        assert result['python_code_sha256']==plan['common_code_sha256']['run_manifest_group_verification.py']
        audit_path=directory/f'd{d}_source_audit.json'
        audit=load(audit_path)
        assert audit['source_sha256']==source['source_sha256']
        assert audit['generator_sha256']==plan['common_code_sha256']['phi_stream_selected.cpp']
        assert audit['candidate_manifest_sha256']==batch['manifest_sha256']
        assert audit['source_row_count']==4*d-2
        assert audit['all_three_moment_relations'] and audit['all_rows_legal_and_complete']
        assert audit['row_order']==source['row_order']
        gp=directory/f'd{d}_gmp_determinants.jsonl'
        xp=directory/f'd{d}_cross_verified_certificates.jsonl'
        gmp,cross=read_lines(gp),read_lines(xp)
        assert cross==result['records']
        assert len(gmp)==len(cross)==result['candidate_count']==batch['candidate_count']
        assert [r['candidate'] for r in gmp]==[r['candidate'] for r in cross]==batch['candidates']
        counts=collections.Counter()
        for i in range(len(cross)):
            r,g=cross[i],gmp[i]
            assert all(r[k]==v for k,v in g.items())
            n=int(r['candidate'])
            assert n in expected and expected[n]==d and n not in actual
            assert r['d']==d and r['N']==4*d and r['source_sha256']==source['source_sha256']
            assert n>8*d and n%(2*d) not in (1,2*d-1)
            assert r['characteristic_gate'] and r['congruence_gate']
            status=r['accepted_status']
            if status=='IF_PRIME_THEN_CONNECTED':
                det=int(r['determinant_mod_n'])
                assert r['status']=='nonzero_determinant' and 0<det<n and math.gcd(det,n)==1
                assert det==int(r['flint_determinant_mod_n'])
                assert r['unit_pivots']==4*d-2 and r['independent_determinants_agree'] is True
            elif status=='EXACT_COMPOSITE_FACTOR':
                f=int(r['factor'])
                assert r['status']=='composite_factor' and 1<f<n and n%f==0
            else:
                assert status=='UNRESOLVED_ZERO_SELECTED_MINOR' and r['status']=='zero_determinant'
                assert int(r['determinant_mod_n'])==int(r['flint_determinant_mod_n'])==0
            counts[status]+=1
            actual[n]=dict(r,complete_receipt_sha256=sha(path))
        assert result['conditional_connected']==counts['IF_PRIME_THEN_CONNECTED']
        assert result['composite']==counts['EXACT_COMPOSITE_FACTOR']
        assert result['unresolved']==counts['UNRESOLVED_ZERO_SELECTED_MINOR']
        receipts.append({str(p.relative_to(ROOT)):sha(p) for p in (path,gp,xp,audit_path,manifest_path)})
    return actual,receipts


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--require-complete',action='store_true');args=ap.parse_args()
    start=time.monotonic()
    accepted_path=ROOT/'control/FINITE_1E56_ROOT_ACCEPTANCE.json'
    accepted=load(accepted_path)
    assert accepted['status']=='ACCEPTED_INTERNAL_PROOF_CHAIN'
    base_path=ROOT/'agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json'
    assert sha(base_path)==accepted['bindings'][str(base_path.relative_to(ROOT))]
    base=load(base_path)
    ds=sorted(base['selected_d']+list(ADDITIONAL))
    assert len(ds)==len(set(ds))==57
    plan,data=load(PLAN),load(DATA)
    assert plan['classification_sha256']==sha(DATA)
    assert plan['campaign_code_sha256']==sha(HERE/'finite_campaign.py')
    assert data['selected_d']==ds and data['lower_exclusive']==str(10**56) and data['upper_inclusive']==str(10**65)
    assert data['inherited_finite_acceptance_sha256']==sha(accepted_path)
    old_code=ROOT/'agents/crt_closure_round3/independent_crt52_closure.py'
    assert sha(old_code)==accepted['bindings'][str(old_code.relative_to(ROOT))]
    spec=importlib.util.spec_from_file_location('audited_other_crt',old_code)
    other=importlib.util.module_from_spec(spec);spec.loader.exec_module(other)
    full=other.enumerate_progressions(ds,10**65)
    values=[n for n in full['candidate_list'] if int(n)>10**56]
    assert values==data['candidate_list']
    assert full['next_candidate_above_upper']==data['next_candidate_above_upper']
    assert full['minimal_modulus']==data['minimal_modulus']
    assert full['sign_components']==data['sign_components']
    assert full['residue_class_count']==data['residue_class_count']
    assert [r['n'] for r in data['records']]==values and len(values)==2689
    pending_general,general_receipts=[],[]
    for d in ADDITIONAL:
        p=ROOT/f'control/D{d}_ROOT_ACCEPTANCE.json'
        if not p.exists():
            pending_general.append(d);continue
        a=load(p)
        assert a['status']=='ACCEPTED_GENERAL_PRIME_D_CRITERION' and a['d']==d and a['N']==4*d
        assert a['condition']==f'p prime, p>{8*d}, p mod {2*d} not in {{1,{2*d-1}}}'
        assert a['all_journal_cells_present'] is True
        for name,digest in a['bindings'].items():
            assert sha(ROOT/name)==digest
        general_receipts.append({'d':d,'acceptance_path':str(p.relative_to(ROOT)),'sha256':sha(p)})
    expected,composite_count={},0
    for r in data['records']:
        n=int(r['n'])
        if r['status']=='EXACT_COMPOSITE_WITNESS':
            check_composite(n,r['certificate']);composite_count+=1
        else:
            assert r['status']=='POINTWISE_EVIDENCE_REQUIRED' and r['primality_proved'] is False
            legal=[d for d in (167,173) if n>8*d and n%(2*d) not in (1,2*d-1)]
            assert legal==r['eligible_d'] and legal[0]==r['preferred_d']
            expected[n]=legal[0]
    assert composite_count==2530 and len(expected)==159
    assert collections.Counter(expected.values())=={167:157,173:2}
    scheduled=[(int(n),b['d']) for b in plan['batches'] for n in b['candidates']]
    assert len(scheduled)==len(set(n for n,d in scheduled))==159 and dict(scheduled)==expected
    sources={d:source_check(d) for d in (167,173)}
    actual,receipts=completed_pointwise(plan,expected,sources)
    missing=sorted(set(expected)-set(actual))
    unresolved=sorted(n for n,r in actual.items() if r['accepted_status']=='UNRESOLVED_ZERO_SELECTED_MINOR')
    ready=not pending_general and not missing and not unresolved
    records=[]
    for r in data['records']:
        n=int(r['n'])
        if r['status']=='EXACT_COMPOSITE_WITNESS':
            records.append(dict(r,verification='Every direct modular-power value recomputed.'))
        elif n not in actual:
            records.append({'n':r['n'],'status':'MISSING_COMPLETED_POINTWISE_RECEIPT','d':expected[n]})
        else:
            a=actual[n]
            records.append({'n':r['n'],'status':a['accepted_status'],'d':expected[n],
                'source_sha256':a['source_sha256'],'determinant_mod_n':a['determinant_mod_n'],
                'complete_receipt_sha256':a['complete_receipt_sha256'],'primality_proved':False})
    out={'status':'ALL_PREREQUISITES_COMPLETE_PENDING_ROOT_ACCEPTANCE' if ready else 'CONDITIONAL_FINITE_CLOSURE_WITH_EXPLICIT_MISSING_GATES',
         'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'interval_accepted':False,
         'closed_interval_certificate_ready':ready,'inherited_upper_inclusive':str(10**56),
         'proposed_upper_inclusive':str(10**65),'candidate_count':2689,'exact_composite_rechecks':composite_count,
         'pointwise_expected':159,'pointwise_completed':len(actual),'pending_general_parameters':pending_general,
         'missing_pointwise_count':len(missing),'unresolved_count':len(unresolved),
         'missing_pointwise':list(map(str,missing)),'unresolved':list(map(str,unresolved)),
         'general_acceptance_receipts':general_receipts,'completed_pointwise_receipts':receipts,
         'source_audits':{str(d):None if a is None else {'sha256':sha(HERE/f'd{d}_INDEPENDENT_SOURCE_AUDIT.json'),'coefficient_count':a['coefficient_count']} for d,a in sources.items()},
         'records':records,'bindings':{str(p.relative_to(ROOT)):sha(p) for p in (accepted_path,base_path,old_code,DATA,PLAN,HERE/'finite_campaign.py')},
         'closure_code_sha256':sha(__file__),'seconds':time.monotonic()-start,
         'scope':'Fresh alternate CRT reconstruction, all composite witnesses, all stored source coefficients/moments and complete dual-backend receipt consistency. No determinant is rerun by this script. Prior small-prime and rank-to-pure-Vieta theorems are inherited. A probable-prime pass is never accepted as primality. Final interval promotion belongs to root.'}
    save(HERE/'FINITE_1E65_CLOSURE_REVIEW.json',out)
    print(json.dumps({k:out[k] for k in ('status','closed_interval_certificate_ready','candidate_count','exact_composite_rechecks','pointwise_completed','pending_general_parameters','missing_pointwise_count','unresolved_count','seconds')},indent=2))
    if args.require_complete:
        assert ready,'General criterion or pointwise evidence is still missing; no interval promotion.'


if __name__=='__main__':main()
