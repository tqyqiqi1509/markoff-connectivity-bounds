#!/usr/bin/env python3
"""Read-only independent finite-chain audit; output only within this package.

No predecessor imports, determinant reruns, claimed primality of candidates,
or changes to root acceptance. Rebuilds a 49,152-entry meet-in-the-middle CRT,
not the full 536,870,912-class Cartesian product.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
from datetime import datetime, timezone
from bisect import bisect_right
import argparse, hashlib, json, math, re, sys, time

sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OLD=ROOT/'agents/finite_extension_round15'
NEW=ROOT/'agents/finite_extension_round20'
MATH={
 'run_manifest_group_verification.py':'286dac5a325ffe9182fade7ff1732b80cdfeded61ddebf8599df493e237853eb',
 'independent_modular_det.cpp':'5d6ea13aaeaf059fd38f7e953e95a105907a733ecbb9e0414199928c8e0e7d31',
 'phi_stream_selected.cpp':'bc0724e42b707c7a269c70b89299013c0cbeb35751f982ba82b9834fb6cdb8d2',
 'phi_stream_partitioned.cpp':'011f27ee78cc2a8ba105429c5b48ab18e58da6a1e538c1e04ff2df87f2e4d851',
 'independent_modular_det':'46fb73c73880e3fd5c31c2d18f354a7bfa39e58720d63ba5ab6e8c2839d5f5c9'}
FIXED={
 'agents/finite_extension_round15/finite_campaign.py':'29fa8018e653badd52d15d1eccc462381ed5d0d699486d19820e186c66f6d49e',
 'agents/finite_extension_round15/BATCH_PLAN.json':'def079ac09e09c1511a1a1cfb4ecf653d7169e85433f7bd2c8dbcfc2c8b9d0f1',
 'agents/finite_extension_round15/CRT57_1E65_INDEPENDENT_CLASSIFICATION.json':'3271f561c259d9fc02d6fc4e66886d12dfe520647e85986a74102201ca25c487',
 'agents/finite_extension_round15/verify_finite_closure.py':'666636c7830930e8ca943937c0afad0ac5c18921e9a68f0448748cc983d7e366',
 'agents/finite_extension_round20/finite_campaign.py':'8c3f1fa967d34e005e691cc06fbf6bc9309fcd06022d0afbae9c24fe7a8857a1',
 'agents/finite_extension_round20/BATCH_PLAN.json':'58431c7806e10574c28cd340d0c69af1b66ce4f1a56aa35bbcb487a3c3ac94a3',
 'agents/finite_extension_round20/CRT57_1E66_INDEPENDENT_CLASSIFICATION.json':'240bd690f918eab6a7dcc50c23ae53df8c822e82ff0033905f23ee35f40386a9',
 'agents/finite_extension_round20/verify_finite_closure.py':'05c1a4cc2ab021b9bf0b81e28452c312ff3000ce11a3afabb7bd7e9d7dfc4c52',
 'agents/parallel_modular_verifier_round40/accept_completed_parallel_uniform.py':'9fec5fcb45e3de7573a530b2ada182c1cba7bc5fd9a759d1665e90af2b13d98a',
 'agents/parallel_modular_verifier_round40/verify_universal_modular_parallel.py':'b8d8a0f34c95df4f8dc1e2a0b87c7ee494050542769f12f41351d53fbc8e2f5e'}
CACHE={};BINDINGS={}

def sha(p):
 p=Path(p);a=p.stat();key=(a.st_dev,a.st_ino,a.st_size,a.st_mtime_ns)
 if key not in CACHE:
  with p.open('rb') as f: digest=hashlib.file_digest(f,'sha256').hexdigest()
  b=p.stat();assert (b.st_dev,b.st_ino,b.st_size,b.st_mtime_ns)==key
  CACHE[key]=digest
 digest=CACHE[key];BINDINGS[str(p.relative_to(ROOT))]={'sha256':digest,'bytes':a.st_size}
 return digest

def load(p):
 p=Path(p);before=sha(p);x=json.loads(p.read_text());assert sha(p)==before;return x

def jsonl(p):
 p=Path(p);before=sha(p);out=[]
 with p.open() as f:
  for line in f:
   assert line.endswith('\n');out.append(json.loads(line))
 assert sha(p)==before;return out

def save(name,obj):
 p=HERE/name;assert not p.exists(),('immutable audit output exists',str(p))
 p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def factor(n):
 ans={};q=2
 while q*q<=n:
  if n%q==0:
   e=0
   while n%q==0:n//=q;e+=1
   ans[q]=e
  q+=1
 if n>1:ans[n]=1
 return ans

def rebuild(ds,lower,upper):
 modulus=math.lcm(*(2*d for d in ds));exps=factor(modulus);parent={q:q for q in exps}
 def find(q):
  while parent[q]!=q:q=parent[q]
  return q
 for d in ds:
  support=[q for q,e in factor(2*d).items() if q!=2 or e>=2]
  for q in support[1:]:parent[find(q)]=find(support[0])
 parts=sorted([sorted(q for q in exps if find(q)==root) for root in set(map(find,exps))])
 blocks=[math.prod(q**exps[q] for q in group) for group in parts]
 assert math.prod(blocks)==modulus and all(b>2 for b in blocks)
 assert all(math.gcd(a,b)==1 for i,a in enumerate(blocks) for b in blocks[i+1:])
 atoms=[]
 for b in blocks:
  m=modulus//b;e=(m*pow(m,-1,b))%modulus
  assert e%b==1 and all(e%c==0 for c in blocks if c!=b)
  atoms.append(e)
 def sums(items):
  out=[0]
  for x in items:out=[(a+x)%modulus for a in out]+[(a-x)%modulus for a in out]
  assert len(out)==len(set(out));return sorted(out)
 k=len(atoms)//2;left=sums(atoms[:k]);right=sums(atoms[k:])
 assert len(left)==16384 and len(right)==32768 and modulus>upper
 vals=[];nexts=[]
 for a in left:
  for shift in (0,modulus):
   lo=bisect_right(right,shift+lower-a);hi=bisect_right(right,shift+upper-a)
   vals.extend(a+b-shift for b in right[lo:hi])
   if hi<len(right):nexts.append(a+right[hi]-shift)
 vals.sort();assert len(vals)==len(set(vals)) and all(lower<n<=upper for n in vals)
 assert all(all(n%(2*d) in (1,2*d-1) for d in ds) for n in vals)
 return {'modulus':modulus,'components':parts,'block_moduli':blocks,
         'class_count':2**len(parts),'half_counts':[len(left),len(right)],
         'values':vals,'next':min(n for n in nexts if n>upper)}

def composite(n,w):
 if w['kind']=='proper_factor':
  a,b=int(w['factor']),int(w['cofactor']);assert 1<a<n and a*b==n;return
 assert w['kind']=='strong_Miller_Rabin_failure'
 a,d,s=int(w['base']),int(w['odd_exponent']),int(w['two_adic_exponent'])
 assert 1<a<n and math.gcd(a,n)==1 and d>0 and d%2==1 and s>0
 assert d*2**s==n-1
 residues=[pow(a,d<<j,n) for j in range(s)]
 assert residues==list(map(int,w['power_chain'])) and residues[0]!=1 and n-1 not in residues

def inspect_source(d):
 r1=OLD/f'd{d}_SOURCE_READY.json';r2=NEW/f'd{d}_SOURCE_READY.json'
 assert sha(r1)==sha(r2);ready=load(r1)
 assert ready['status']=='COMPLETE_SOURCE_IDENTITY_BOUND_NOT_A_RANK_CLAIM'
 assert ready['d']==d and ready['N']==4*d
 assert ready['generator_sha256']==MATH['phi_stream_partitioned.cpp']
 assert ready['legacy_selected_generator_equivalent_prefix_sha256']==MATH['phi_stream_selected.cpp']
 for name,s in ready['bindings'].items():assert sha(ROOT/name)==s
 src=ROOT/ready['source_path'];assert sha(src)==ready['source_sha256'] and src.stat().st_size==ready['source_bytes']
 for name in ('phi_stream_selected.cpp','phi_stream_partitioned.cpp'):
  assert sha(src.parent/name)==MATH[name]
 expected=[(j*d,m) for j in range(1,5) for m in range(j*(d-1)//2 if j<4 else d+1)]
 assert len(expected)==4*d-2
 metadata={};count=0;maxdegree=0;binom=[math.comb(2*j,j) for j in range(4*d+1)]
 with src.open() as f:
  for line_no,line in enumerate(f):
   assert line.endswith('\n');row=json.loads(line);key=(row['n'],row['m'])
   assert key in expected and key not in metadata and row['d']==d and row['nd']==4*d
   strings=row['coefficients'];values=list(map(int,strings))
   assert len(values)==key[0]+1 and all(str(v)==s for v,s in zip(values,strings))
   assert values[0]==0 and values[1]==12*values[2]
   assert sum(v*binom[j] for j,v in enumerate(values))==0
   degree=2*key[0]+(key[0]//d)//2+key[1];assert row['max_half_total_degree']==degree
   maxdegree=max(maxdegree,degree);count+=len(values)
   metadata[key]={'matrix_row':expected.index(key),'source_line_zero_based':line_no,'n':key[0],'m':key[1],
     'full_coefficients_sha256':hashlib.sha256((','.join(strings)+'\n').encode()).hexdigest()}
 assert set(metadata)==set(expected) and sha(src)==ready['source_sha256']
 out={'d':d,'N':4*d,'source_sha256':ready['source_sha256'],'source_bytes':ready['source_bytes'],
      'source_ready_sha256':sha(r1),'rows':len(expected),'coefficient_count':count,
      'max_half_total_degree':maxdegree,'row_order':[metadata[k] for k in expected],
      'scope':'Fresh full integer parsing, all three moment identities and row keys; no independent polynomial regeneration or determinant computation.'}
 return out

def pointwise(directory,plan,class_path,expected,sources):
 records={};missing=[];summary=[]
 assert plan['classification_sha256']==sha(class_path)
 assert plan['campaign_code_sha256']==sha(directory/'finite_campaign.py')
 assert plan['common_code_sha256']==MATH
 schedule=[(int(n),b['d']) for b in plan['batches'] for n in b['candidates']]
 assert len(schedule)==len(set(n for n,d in schedule))==len(expected) and dict(schedule)==expected
 for b in plan['batches']:
  folder=directory/b['directory'];d=b['d'];src=sources[d];manifest_path=folder/'batch_manifest.json'
  assert sha(manifest_path)==b['manifest_sha256'];manifest=load(manifest_path)
  assert manifest['groups']==[{'d':d,'N':4*d,'candidate_count':b['candidate_count'],'candidates':b['candidates']}]
  assert manifest['classification_sha256']==sha(class_path)
  assert (folder/f'd{d}_candidates.txt').read_text().splitlines()==b['candidates']
  sha(folder/f'd{d}_candidates.txt')
  for name,digest in MATH.items():assert sha(folder/name)==digest
  assert (folder/'run_manifest_group_verification.py').resolve().parent==folder.resolve()
  assert sha(folder/f'd{d}_prefix_rows.jsonl')==src['source_sha256']
  final_path=folder/f'd{d}_complete_verification.json'
  if not final_path.exists():missing.extend(b['candidates']);continue
  final=load(final_path);audit=load(folder/f'd{d}_source_audit.json')
  assert final['status']=='COMPLETE' and final['d']==d and final['N']==4*d
  assert final['source_sha256']==src['source_sha256']
  assert final['gmp_code_sha256']==MATH['independent_modular_det.cpp']
  assert final['python_code_sha256']==MATH['run_manifest_group_verification.py']
  assert audit['source_sha256']==src['source_sha256'] and audit['source_bytes']==src['source_bytes']
  assert audit['generator_sha256']==MATH['phi_stream_selected.cpp']
  assert audit['candidate_manifest_sha256']==sha(manifest_path)
  assert audit['source_row_count']==src['rows'] and audit['all_three_moment_relations'] is True and audit['all_rows_legal_and_complete'] is True
  assert audit['row_order']==src['row_order']
  gmp=jsonl(folder/f'd{d}_gmp_determinants.jsonl');cross=jsonl(folder/f'd{d}_cross_verified_certificates.jsonl')
  assert cross==final['records']
  assert len(gmp)==len(cross)==final['candidate_count']==b['candidate_count']
  assert [r['candidate'] for r in gmp]==[r['candidate'] for r in cross]==b['candidates']
  counts=Counter()
  for g,c in zip(gmp,cross):
   assert all(c[k]==v for k,v in g.items());n=int(c['candidate'])
   assert n not in records and expected[n]==d and c['d']==d and c['N']==4*d
   assert c['source_sha256']==src['source_sha256']
   assert n>max(8*d,src['max_half_total_degree']) and n%(2*d) not in (1,2*d-1)
   assert c['characteristic_gate'] is True and c['congruence_gate'] is True
   status=c['accepted_status'];counts[status]+=1
   if status=='IF_PRIME_THEN_CONNECTED':
    v=int(c['determinant_mod_n']);assert c['status']=='nonzero_determinant' and 0<v<n
    assert math.gcd(v,n)==1 and v==int(c['flint_determinant_mod_n'])
    assert c['unit_pivots']==4*d-2 and c['independent_determinants_agree'] is True
   elif status=='EXACT_COMPOSITE_FACTOR':
    f=int(c['factor']);assert c['status']=='composite_factor' and 1<f<n and n%f==0
   else:
    assert status=='UNRESOLVED_ZERO_SELECTED_MINOR' and c['status']=='zero_determinant'
    assert int(c['determinant_mod_n'])==int(c['flint_determinant_mod_n'])==0
   records[n]={'candidate':str(n),'d':d,'status':status,'determinant_mod_n':c.get('determinant_mod_n'),
    'source_sha256':src['source_sha256'],'complete_receipt_sha256':sha(final_path),
    'batch':str(folder.relative_to(ROOT)),'primality_proved':False}
  assert final['conditional_connected']==counts['IF_PRIME_THEN_CONNECTED']
  assert final['composite']==counts['EXACT_COMPOSITE_FACTOR'] and final['unresolved']==counts['UNRESOLVED_ZERO_SELECTED_MINOR']
  summary.append({'batch':str(folder.relative_to(ROOT)),'d':d,'count':len(cross),'counts':dict(counts),'receipt_sha256':sha(final_path)})
 return records,missing,summary

def general_snapshot():
 ready=[];pending=[]
 for d in [139,149,151,157,163]:
  p=ROOT/f'control/D{d}_ROOT_ACCEPTANCE.json'
  if not p.exists():pending.append(d);continue
  a=load(p)
  assert a['status']=='ACCEPTED_GENERAL_PRIME_D_CRITERION' and a['d']==d and a['N']==4*d
  assert a['condition']==f'p prime, p>{8*d}, p mod {2*d} not in {{1,{2*d-1}}}'
  assert a['all_journal_cells_present'] is True
  for name,digest in a['bindings'].items():assert sha(ROOT/name)==digest
  if 'parallel_implementation_id' in a:
   assert a['parallel_implementation_id']=='proven_proth_parallel_nmod_v1'
   assert a['root_closeout_code_sha256']==FIXED['agents/parallel_modular_verifier_round40/accept_completed_parallel_uniform.py']
  ready.append({'d':d,'path':str(p.relative_to(ROOT)),'sha256':sha(p),'binding_count':len(a['bindings'])})
 return {'accepted':ready,'pending':pending}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--label',required=True);args=ap.parse_args()
 assert re.fullmatch('[a-zA-Z0-9_]+',args.label);started=time.monotonic()
 for name,digest in FIXED.items():assert sha(ROOT/name)==digest
 old_data=load(OLD/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json');data=load(NEW/'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json')
 old_plan=load(OLD/'BATCH_PLAN.json');plan=load(NEW/'BATCH_PLAN.json')
 accepted_path=ROOT/'control/FINITE_1E56_ROOT_ACCEPTANCE.json';acc=load(accepted_path)
 assert acc['status']=='ACCEPTED_INTERNAL_PROOF_CHAIN'
 for name,digest in acc['bindings'].items():assert sha(ROOT/name)==digest
 base=load(ROOT/'agents/crt_closure_round3/CRT52_INDEPENDENT_ENUMERATION.json')
 ds=sorted(base['selected_d']+[139,149,151,157,163]);assert ds==data['selected_d'] and len(ds)==len(set(ds))==57
 assert data['inherited_finite_acceptance_sha256']==sha(accepted_path)
 assert data['lower_exclusive']==str(10**56) and data['upper_inclusive']==str(10**66)
 crt=rebuild(ds,10**56,10**66);vals=crt.pop('values')
 assert list(map(str,vals))==data['candidate_list'] and len(vals)==27985
 assert str(crt['modulus'])==data['minimal_modulus'] and crt['components']==data['sign_components']
 assert crt['class_count']==data['residue_class_count'] and str(crt['next'])==data['next_candidate_above_upper']
 sequence_sha=hashlib.sha256(('\n'.join(map(str,vals))+'\n').encode()).hexdigest()
 assert sequence_sha==data['candidate_sequence_sha256']=='5399b99b3aace9bf17fcf322d40dc8eaf3e2b8a2cb2cd7c9eb5a3c5afa4c9809'
 assert [r['n'] for r in data['records']]==data['candidate_list']
 assert [r for r in data['records'] if int(r['n'])<=10**65]==old_data['records']
 assert [str(n) for n in vals if n<=10**65]==old_data['candidate_list']
 assert str(next(n for n in vals if n>10**65))==old_data['next_candidate_above_upper']
 expected={};composites=0
 for r in data['records']:
  n=int(r['n'])
  if r['status']=='EXACT_COMPOSITE_WITNESS':composite(n,r['certificate']);composites+=1
  else:
   assert r['status']=='POINTWISE_EVIDENCE_REQUIRED' and r['primality_proved'] is False
   eligible=[d for d in [167,173] if n>8*d and n%(2*d) not in (1,2*d-1)]
   assert eligible==r['eligible_d'] and eligible[0]==r['preferred_d'];expected[n]=eligible[0]
 assert composites==26210 and Counter(expected.values())=={167:1755,173:20}
 assert list(map(str,expected))==data['pointwise_candidates']
 reused={int(n):expected[int(n)] for n in old_data['pointwise_candidates']}
 fresh={n:d for n,d in expected.items() if n not in reused}
 assert len(reused)==159 and Counter(reused.values())=={167:157,173:2}
 assert len(fresh)==1616 and Counter(fresh.values())=={167:1598,173:18}
 assert plan['reuse_dependency']['candidates']==old_data['pointwise_candidates']
 for name,digest in plan['reuse_dependency']['static_bindings'].items():assert sha(ROOT/name)==digest
 wl_path=NEW/'WOLFRAM_CRT57_1E66.wl';wl=wl_path.read_text();sha(wl_path)
 assert list(map(int,re.search(r'ds=\{([^}]+)\}',wl).group(1).split(',')))==ds
 assert list(map(int,re.search(r'targets=\{([^}]+)\}',wl).group(1).split(',')))==[n for n,d in expected.items() if d==173]
 comp=load(NEW/'WOLFRAM_CRT57_1E66_COMPARISON.json')
 for name,digest in comp['bindings'].items():assert sha(NEW/name)==digest
 raw=load(NEW/'WOLFRAM_CRT57_1E66_RAW.json');assert raw['isError'] is False
 text='\n'.join(c['text'] for c in raw['content'] if c['type']=='text')
 for key,value in {'modulus':crt['modulus'],'signComponents':29,'classes':2**29,'count':27985,'sum':sum(vals),'sumSquares':sum(n*n for n in vals)}.items():
  assert int(re.search('"'+key+r'" -> ([0-9]+)',text).group(1))==value
 assert re.search(r'"sequenceSHA256" -> "([0-9a-f]+)"',text).group(1)==sequence_sha
 assert re.search(r'"hashABC" -> "([0-9a-f]+)"',text).group(1)==hashlib.sha256(b'abc').hexdigest()
 for key in ['unique','allCongruences','source173Only']:assert '"'+key+'" -> True' in text
 sources={d:inspect_source(d) for d in [167,173]}
 old_actual,old_missing,old_batches=pointwise(OLD,old_plan,OLD/'CRT57_1E65_INDEPENDENT_CLASSIFICATION.json',reused,sources)
 new_actual,new_missing,new_batches=pointwise(NEW,plan,NEW/'CRT57_1E66_INDEPENDENT_CLASSIFICATION.json',fresh,sources)
 assert not set(old_actual)&set(new_actual);actual={**old_actual,**new_actual}
 missing=sorted(set(expected)-set(actual));assert missing==sorted(map(int,old_missing+new_missing))
 unresolved=[n for n,r in actual.items() if r['status']=='UNRESOLVED_ZERO_SELECTED_MINOR']
 general=general_snapshot();ready=not missing and not unresolved and not general['pending']
 out={'status':'ALL_FINITE_DEPENDENCIES_VERIFIED_PENDING_ROOT_ACCEPTANCE' if ready else 'CONDITIONAL_REVIEW_WITH_EXPLICIT_MISSING_GATES',
      'utc':datetime.now(timezone.utc).isoformat(),'interval_accepted':False,'proposed_upper_inclusive':str(10**66),
      'inherited_upper_inclusive':str(10**56),'crt':crt,'candidate_count':len(vals),'sequence_sha256':sequence_sha,
      'fresh_composite_checks':composites,'pointwise_expected':1775,'pointwise_complete':len(actual),
      'reused_complete':len(old_actual),'fresh_complete':len(new_actual),'completed_batches':len(old_batches)+len(new_batches),
      'missing_count':len(missing),'missing':list(map(str,missing)),'unresolved_count':len(unresolved),
      'general':general,'full_interval_certificate_ready':ready,
      'intermediate_1e65':{'candidate_count':2689,'pointwise_expected':159,'pointwise_complete':len(old_actual),
       'missing_count':len(old_missing),'unresolved_count':sum(n<=10**65 for n in unresolved),
       'ready':not old_missing and not any(n<=10**65 for n in unresolved) and not general['pending']},
      'pointwise_status_counts':dict(Counter(x['status'] for x in actual.values())),
      'sources':sources,'batch_summaries':old_batches+new_batches,'records':[actual[n] for n in sorted(actual)],
      'input_bindings':BINDINGS,'reviewer_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'seconds':time.monotonic()-started,
      'scope':'No predecessor imports or mutations. Fresh two-half CRT reconstruction, all direct modular composite powers, all source integer moments and all available atomic dual-backend receipts. No new determinant/poly generation; prior graph theorem and root general acceptances remain explicit prerequisites.'}
 save('REVIEW_'+args.label+'.json',out)
 print(json.dumps({k:v for k,v in out.items() if k not in ['input_bindings','sources','batch_summaries','records','crt']},ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
