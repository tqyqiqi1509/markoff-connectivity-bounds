#!/usr/bin/env python3
"""Independent arithmetic/source/receipt closure, without new determinants.

Rebuilds CRT using the older half-modulus progression implementation, distinct
from this branch's full-modulus signed-atom producer. Every composite power is
recomputed directly rather than by the producer's iterative squaring chain.
All source coefficients and moments are parsed/rechecked once per ready d.
Root remains responsible for final mathematical acceptance.
"""
from pathlib import Path
import argparse, collections, datetime, hashlib, importlib.util, json, math, time, sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE/'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json'
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
        ready_path=ROOT/'agents/finite_extension_round15'/f'd{d}_SOURCE_READY.json'
    if not ready_path.exists():
        return None
    ready=load(ready_path)
    assert ready['d']==d and ready['N']==4*d
    assert ready['actual_generator']=='phi_stream_partitioned.cpp'
    assert ready['generator_sha256']=='011f27ee78cc2a8ba105429c5b48ab18e58da6a1e538c1e04ff2df87f2e4d851'
    assert ready['legacy_selected_generator_equivalent_prefix_sha256']=='bc0724e42b707c7a269c70b89299013c0cbeb35751f982ba82b9834fb6cdb8d2'
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
            'source_ready_path':str(ready_path.relative_to(ROOT)), 'source_ready_sha256':sha(ready_path),'scope':'Fresh complete numeric parsing and moment/row-order checks; not a second polynomial generation or determinant run.'}
    save(HERE/f'd{d}_INDEPENDENT_SOURCE_AUDIT.json',result)
    return result


def read_lines(path):
    result=[]
    with path.open() as f:
        for line in f:
            assert line.endswith('\n')
            result.append(json.loads(line))
    return result


def completed_pointwise(plan,expected,sources,campaign_dir,classification_path):
    actual,receipts={},[]
    for batch in plan['batches']:
        d=batch['d'];directory=campaign_dir/batch['directory']
        path=directory/f'd{d}_complete_verification.json'
        if not path.exists():
            continue
        assert sources[d] is not None
        source=sources[d]
        ready_path=campaign_dir/f'd{d}_SOURCE_READY.json'
        assert sha(ready_path)==source['source_ready_sha256']
        manifest_path=directory/'batch_manifest.json'
        assert sha(manifest_path)==batch['manifest_sha256']
        manifest=load(manifest_path)
        assert len(manifest['groups'])==1
        group=manifest['groups'][0]
        assert group['d']==d and group['N']==4*d
        assert group['candidates']==batch['candidates']
        assert group['candidate_count']==len(group['candidates'])==batch['candidate_count']
        assert manifest['classification_sha256']==sha(classification_path)
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
        receipts.append({str(p.relative_to(ROOT)):sha(p) for p in (path,gp,xp,audit_path,manifest_path,ready_path)})
    return actual,receipts


REUSE_DIR=ROOT/'agents/finite_extension_round15'
REUSE_FIXED={
    'finite_campaign.py':'29fa8018e653badd52d15d1eccc462381ed5d0d699486d19820e186c66f6d49e',
    'BATCH_PLAN.json':'def079ac09e09c1511a1a1cfb4ecf653d7169e85433f7bd2c8dbcfc2c8b9d0f1',
    'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json':'3271f561c259d9fc02d6fc4e66886d12dfe520647e85986a74102201ca25c487',
    'verify_finite_closure.py':'666636c7830930e8ca943937c0afad0ac5c18921e9a68f0448748cc983d7e366',
}
MATHEMATICS={
    'run_manifest_group_verification.py':'286dac5a325ffe9182fade7ff1732b80cdfeded61ddebf8599df493e237853eb',
    'independent_modular_det.cpp':'5d6ea13aaeaf059fd38f7e953e95a105907a733ecbb9e0414199928c8e0e7d31',
    'phi_stream_selected.cpp':'bc0724e42b707c7a269c70b89299013c0cbeb35751f982ba82b9834fb6cdb8d2',
    'phi_stream_partitioned.cpp':'011f27ee78cc2a8ba105429c5b48ab18e58da6a1e538c1e04ff2df87f2e4d851',
    'independent_modular_det':'46fb73c73880e3fd5c31c2d18f354a7bfa39e58720d63ba5ab6e8c2839d5f5c9',
}


def inherited_reuse(plan,data,expected):
    required={str((REUSE_DIR/name).relative_to(ROOT)):digest for name,digest in REUSE_FIXED.items()}
    assert plan['reuse_dependency']['static_bindings']==required
    for name,digest in required.items():
        assert sha(ROOT/name)==digest
    old_data_path=REUSE_DIR/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json'
    old_data=load(old_data_path);old_plan=load(REUSE_DIR/'BATCH_PLAN.json')
    assert old_plan['classification_sha256']==sha(old_data_path)
    assert old_plan['campaign_code_sha256']==sha(REUSE_DIR/'finite_campaign.py')
    assert old_plan['common_code_sha256']==plan['common_code_sha256']==MATHEMATICS
    assert [n for n in data['candidate_list'] if int(n)<=10**65]==old_data['candidate_list']
    assert [r for r in data['records'] if int(r['n'])<=10**65]==old_data['records']
    assert old_data['candidate_count']==2689 and old_data['pointwise_count']==159
    assert plan['reuse_dependency']['candidates']==old_data['pointwise_candidates']
    reused={int(n):expected[int(n)] for n in old_data['pointwise_candidates']}
    assert len(reused)==159 and collections.Counter(reused.values())=={167:157,173:2}
    schedule=[(int(n),b['d']) for b in old_plan['batches'] for n in b['candidates']]
    assert len(schedule)==len(set(n for n,d in schedule))==159 and dict(schedule)==reused
    for batch in old_plan['batches']:
        assert sha(REUSE_DIR/batch['directory']/'batch_manifest.json')==batch['manifest_sha256']
    return old_plan,old_data_path,reused,required


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--require-complete',action='store_true');args=ap.parse_args()
    start=time.monotonic()
    accepted_path=ROOT/'control/FINITE_1E56_ROOT_ACCEPTANCE.json';accepted=load(accepted_path)
    assert accepted['status']=='ACCEPTED_INTERNAL_PROOF_CHAIN'
    for name,digest in accepted['bindings'].items():
        assert sha(ROOT/name)==digest
    inherited_closure=load(ROOT/'agents/crt_closure_round3/CRT52_FINITE_CLOSURE_REVIEW.json')
    assert inherited_closure['interval_accepted'] is True
    assert inherited_closure['upper_inclusive']==str(10**56)
    assert inherited_closure['missing_count']==inherited_closure['unresolved_count']==0
    base_path=ROOT/'agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json';base=load(base_path)
    ds=sorted(base['selected_d']+list(ADDITIONAL));assert len(ds)==len(set(ds))==57
    plan,data=load(PLAN),load(DATA)
    assert plan['classification_sha256']==sha(DATA)
    assert plan['campaign_code_sha256']==sha(HERE/'finite_campaign.py')
    assert plan['common_code_sha256']==MATHEMATICS
    binary_binding=load(ROOT/'control/FINITE_BINARY_ROOT_BINDING.json')
    assert binary_binding['binary_sha256']==MATHEMATICS['independent_modular_det']
    assert data['selected_d']==ds and data['lower_exclusive']==str(10**56) and data['upper_inclusive']==str(10**66)
    assert data['inherited_finite_acceptance_sha256']==sha(accepted_path)
    assert data['code_sha256']==sha(HERE/'rebuild_crt57.py')
    old_code=ROOT/'agents/crt_closure_round3/independent_crt52_closure.py'
    assert sha(old_code)==accepted['bindings'][str(old_code.relative_to(ROOT))]
    spec=importlib.util.spec_from_file_location('audited_other_crt',old_code)
    other=importlib.util.module_from_spec(spec);spec.loader.exec_module(other)
    full=other.enumerate_progressions(ds,10**66)
    values=[n for n in full['candidate_list'] if int(n)>10**56]
    assert values==data['candidate_list']
    for key in ('next_candidate_above_upper','minimal_modulus','sign_components','residue_class_count'):
        assert full[key]==data[key],key
    sequence_sha=hashlib.sha256(('\n'.join(values)+'\n').encode()).hexdigest()
    assert sequence_sha==data['candidate_sequence_sha256']
    assert [r['n'] for r in data['records']]==values and len(values)==27985
    # The 10^65 intermediate boundary is independently reconstructed too.
    at65=other.enumerate_progressions(ds,10**65)
    old_class=load(REUSE_DIR/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json')
    assert [n for n in at65['candidate_list'] if int(n)>10**56]==old_class['candidate_list']
    assert at65['next_candidate_above_upper']==old_class['next_candidate_above_upper']
    pending_general,general_receipts=[],[]
    for d in ADDITIONAL:
        p=ROOT/f'control/D{d}_ROOT_ACCEPTANCE.json'
        if not p.exists():pending_general.append(d);continue
        a=load(p)
        assert a['status']=='ACCEPTED_GENERAL_PRIME_D_CRITERION' and a['d']==d and a['N']==4*d
        assert a['condition']==f'p prime, p>{8*d}, p mod {2*d} not in {{1,{2*d-1}}}'
        assert a['all_journal_cells_present'] is True
        for name,digest in a['bindings'].items():assert sha(ROOT/name)==digest
        general_receipts.append({'d':d,'acceptance_path':str(p.relative_to(ROOT)),'sha256':sha(p)})
    expected,composite_count={},0
    for r in data['records']:
        n=int(r['n'])
        if r['status']=='EXACT_COMPOSITE_WITNESS':
            check_composite(n,r['certificate']);composite_count+=1
        else:
            assert r['status']=='POINTWISE_EVIDENCE_REQUIRED' and r['primality_proved'] is False
            legal=[d for d in (167,173) if n>8*d and n%(2*d) not in (1,2*d-1)]
            assert legal==r['eligible_d'] and legal[0]==r['preferred_d'];expected[n]=legal[0]
    assert composite_count==26210 and len(expected)==1775
    assert collections.Counter(expected.values())=={167:1755,173:20}
    old_plan,old_data_path,reused,reuse_bindings=inherited_reuse(plan,data,expected)
    fresh={n:d for n,d in expected.items() if n not in reused}
    assert len(fresh)==1616 and collections.Counter(fresh.values())=={167:1598,173:18}
    scheduled=[(int(n),b['d']) for b in plan['batches'] for n in b['candidates']]
    assert len(scheduled)==len(set(n for n,d in scheduled))==1616 and dict(scheduled)==fresh
    assert plan['planned_total']==1616 and plan['total_pointwise_including_reuse']==1775
    assert plan['workers']==4 and plan['max_batch_size']==20 and len(plan['batches'])==81
    assert all(0<b['candidate_count']<=20 and b['worker']==b['id']%4 for b in plan['batches'])
    assert [b['id'] for b in plan['batches']]==list(range(81))
    sources={d:source_check(d) for d in (167,173)}
    actual,receipts=completed_pointwise(plan,fresh,sources,HERE,DATA)
    old_actual,old_receipts=completed_pointwise(old_plan,reused,sources,REUSE_DIR,old_data_path)
    assert not(set(actual)&set(old_actual));actual.update(old_actual);receipts.extend(old_receipts)
    missing=sorted(set(expected)-set(actual))
    unresolved=sorted(n for n,r in actual.items() if r['accepted_status']=='UNRESOLVED_ZERO_SELECTED_MINOR')
    ready=not pending_general and not missing and not unresolved
    records=[]
    for r in data['records']:
        n=int(r['n'])
        if r['status']=='EXACT_COMPOSITE_WITNESS':
            records.append(dict(r,verification='Every direct modular-power value recomputed.'))
        elif n not in actual:
            records.append({'n':r['n'],'status':'MISSING_COMPLETED_POINTWISE_RECEIPT','d':expected[n],
                            'origin':'round15_pending_reuse' if n in reused else 'round20_new'})
        else:
            a=actual[n]
            record={'n':r['n'],'status':a['accepted_status'],'d':expected[n],
                    'origin':'round15_verified_reuse' if n in reused else 'round20_new',
                    'source_sha256':a['source_sha256'],'complete_receipt_sha256':a['complete_receipt_sha256'],'primality_proved':False}
            for key in ('determinant_mod_n','flint_determinant_mod_n','factor'):
                if key in a:record[key]=a[key]
            records.append(record)
    bad=sorted(set(missing)|set(unresolved))
    prefix_upper=min(10**66,bad[0]-1) if bad else 10**66
    prefix_values=[int(n) for n in values if int(n)<=prefix_upper]
    next_after_prefix=next((n for n in values if int(n)>prefix_upper),data['next_candidate_above_upper'])
    prefix={'upper_inclusive':str(prefix_upper),'first_missing_or_unresolved_pointwise':str(bad[0]) if bad else None,
            'last_necessary_integer':str(prefix_values[-1]) if prefix_values else None,
            'next_necessary_integer':next_after_prefix,
            'arithmetic_and_completed_batch_pointwise_coverage':True,
            'all_general_gates_complete':not pending_general,
            'certificate_ready_pending_root_acceptance':not pending_general,
            'interval_accepted':False,
            'rule':'Only atomic COMPLETE batches enter this prefix; every lesser necessary integer was either strictly composite or has complete matching GMP/FLINT evidence. Pending general gates still forbid extension of the accepted interval.'}
    at65_expected={n:d for n,d in expected.items() if n<=10**65}
    at65_missing=[str(n) for n in missing if n<=10**65];at65_unresolved=[str(n) for n in unresolved if n<=10**65]
    at65_ready=not pending_general and not at65_missing and not at65_unresolved
    intermediate={'status':'ALL_1E65_PREREQUISITES_COMPLETE_PENDING_ROOT_ACCEPTANCE' if at65_ready else 'CONDITIONAL_1E65_INTERMEDIATE_WITH_EXPLICIT_MISSING_GATES',
        'interval_accepted':False,'closed_interval_certificate_ready':at65_ready,
        'lower_exclusive':str(10**56),'proposed_upper_inclusive':str(10**65),
        'last_candidate':old_class['candidate_list'][-1],'next_candidate_above_upper':at65['next_candidate_above_upper'],
        'candidate_count':2689,'pointwise_expected':len(at65_expected),'pointwise_completed':sum(n in actual for n in at65_expected),
        'pending_general_parameters':pending_general,'missing_pointwise':at65_missing,'unresolved':at65_unresolved,
        'inherited_finite_acceptance_sha256':sha(accepted_path),'old_classification_sha256':sha(old_data_path),
        'records':[r for r in records if int(r['n'])<=10**65],'closure_code_sha256':sha(__file__)}
    save(HERE/'FINITE_1E65_INTERMEDIATE_CLOSURE_REVIEW.json',intermediate)
    out={'status':'ALL_PREREQUISITES_COMPLETE_PENDING_ROOT_ACCEPTANCE' if ready else 'CONDITIONAL_FINITE_CLOSURE_WITH_EXPLICIT_MISSING_GATES',
         'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'interval_accepted':False,
         'closed_interval_certificate_ready':ready,'inherited_upper_inclusive':str(10**56),
         'proposed_upper_inclusive':str(10**66),'candidate_count':27985,'exact_composite_rechecks':composite_count,
         'first_candidate_above_inherited_upper':values[0],'last_candidate_at_or_below_upper':values[-1],
         'next_candidate_above_upper':data['next_candidate_above_upper'],'candidate_sequence_sha256':sequence_sha,
         'pointwise_expected':1775,'pointwise_completed':len(actual),'pending_general_parameters':pending_general,
         'reused_pointwise_expected':159,'reused_pointwise_completed':len(old_actual),
         'new_pointwise_expected':1616,'new_pointwise_completed':len(actual)-len(old_actual),
         'old_173_pending_reuse':[str(n) for n in reused if reused[n]==173 and n not in old_actual],
         'missing_pointwise_count':len(missing),'unresolved_count':len(unresolved),
         'missing_pointwise':list(map(str,missing)),'unresolved':list(map(str,unresolved)),
         'maximum_complete_batch_prefix':prefix,'intermediate_1e65':{'sha256':sha(HERE/'FINITE_1E65_INTERMEDIATE_CLOSURE_REVIEW.json'),'ready':at65_ready},
         'general_acceptance_receipts':general_receipts,'completed_pointwise_receipts':receipts,
         'source_audits':{str(d):None if a is None else {'sha256':sha(HERE/f'd{d}_INDEPENDENT_SOURCE_AUDIT.json'),'coefficient_count':a['coefficient_count']} for d,a in sources.items()},
         'records':records,'bindings':{str(p.relative_to(ROOT)):sha(p) for p in (accepted_path,base_path,old_code,DATA,PLAN,HERE/'finite_campaign.py',HERE/'rebuild_crt57.py')},
         'inherited_acceptance_bindings_all_rechecked':True,'reuse_static_bindings':reuse_bindings,
         'closure_code_sha256':sha(__file__),'seconds':time.monotonic()-start,
         'scope':'Fresh alternate CRT reconstruction through both 10^65 and 10^66, every strict composite witness, every stored ready source coefficient/moment, and all available complete dual-backend receipts including 159 explicitly tracked legacy candidates. No determinants are rerun. Source generator identities and inherited 10^56 chain are bound. No probable-prime pass is primality. Root alone can accept an interval.'}
    save(HERE/'FINITE_1E66_CLOSURE_REVIEW.json',out)
    print(json.dumps({k:out[k] for k in ('status','closed_interval_certificate_ready','candidate_count','exact_composite_rechecks','pointwise_completed','reused_pointwise_completed','new_pointwise_completed','pending_general_parameters','missing_pointwise_count','unresolved_count','seconds')},indent=2))
    if args.require_complete:assert ready,'General criterion or pointwise evidence is missing; no interval promotion.'


if __name__=='__main__':main()
