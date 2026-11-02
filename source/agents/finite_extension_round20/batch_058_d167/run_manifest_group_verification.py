#!/usr/bin/env python3
"""Complete integer source checks and two independent modular determinants.
The conditional mathematical conclusion is only: if candidate is prime, its
pure-Vieta Markoff graph is connected. No probabilistic primality claim.
"""
import argparse,datetime,hashlib,json,math,os,subprocess,sys,time
from pathlib import Path
from flint import fmpz_mod_ctx,fmpz_mod_mat
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def save(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');temp.replace(path)

def verify_group(d, group_manifest):
 N=4*d;expected=[(b*d,m) for b in range(1,5) for m in range(b*(d-1)//2 if b<4 else d+1)]
 assert len(expected)==N-2 and all(d%q for q in range(2,math.isqrt(d)+1))
 source=BASE/f'd{d}_prefix_rows.jsonl';log=BASE/f'd{d}_generation.log'
 assert source.is_file(), 'Run the source generator to completion first'
 records={};bins=[math.comb(2*i,i) for i in range(N+1)];line_indices={};canonical=[]
 with source.open() as f:
  for lineno,line in enumerate(f):
   assert line.endswith('\n');r=json.loads(line)
   assert r['d']==d and r['nd']==N
   key=(r['n'],r['m']);assert key in expected and key not in records
   co=list(map(int,r['coefficients']));assert len(co)==r['n']+1
   assert co[0]==0 and co[1]==12*co[2]
   assert sum(c*b for c,b in zip(co,bins))==0
   assert r['max_half_total_degree']==2*r['n']+(r['n']//d)//2+r['m']
   records[key]=co;line_indices[key]=lineno
 assert set(records)==set(expected)
 rows=[];metadata=[]
 for i,key in enumerate(expected):
  co=records[key];rows.append(co[3:]+[0]*(N+1-len(co)))
  assert len(rows[-1])==N-2
  h=hashlib.sha256((','.join(map(str,co))+'\n').encode()).hexdigest()
  metadata.append({'matrix_row':i,'source_line_zero_based':line_indices[key],'n':key[0],'m':key[1],'full_coefficients_sha256':h})
 del records
 grouping=json.load(group_manifest.open())
 matches=[g for g in grouping['groups'] if g['d']==d]
 assert len(matches)==1 and matches[0]['N']==N
 candidates=list(map(int,matches[0]['candidates']))
 assert len(candidates)==len(set(candidates))==matches[0]['candidate_count']
 assert all(p>2*N and p%(2*d) not in (1,2*d-1) for p in candidates)
 assert candidates==list(map(int,(BASE/f'd{d}_candidates.txt').read_text().split()))
 source_audit={'status':'ALL_INTEGER_SOURCE_CHECKS_PASSED','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'d':d,'N':N,'source_file':str(source),'source_sha256':sha(source),'source_bytes':source.stat().st_size,'generator_sha256':sha(BASE/'phi_stream_selected.cpp'),'source_row_count':len(rows),'columns':'coefficients a_3 through a_N for x^(2j), zero padded','all_three_moment_relations':True,'all_rows_legal_and_complete':True,'row_order':metadata,'candidate_count':len(candidates),'scope':'All rows and integer identities checked; no claim of a second independent full polynomial generation.'}
 source_audit['candidate_manifest_file']=str(group_manifest)
 source_audit['candidate_manifest_sha256']=sha(group_manifest)
 save(BASE/f'd{d}_source_audit.json',source_audit)
 print(json.dumps({'stage':'source_audit_passed','d':d,'rows':len(rows),'bytes':source.stat().st_size}),flush=True)
 cppout=BASE/f'd{d}_gmp_determinants.jsonl';cpplog=BASE/f'd{d}_gmp_determinants.log'
 started=time.monotonic()
 with cpplog.open('w') as f:
  run=subprocess.run([str(BASE/'independent_modular_det'),str(d),str(source),str(BASE/f'd{d}_candidates.txt'),str(cppout)],stderr=f)
 assert run.returncode==0,('GMP verifier failed',run.returncode)
 cpp=[json.loads(x) for x in cppout.read_text().splitlines()]
 assert len(cpp)==len(candidates)
 results=[]
 with (BASE/f'd{d}_cross_verified_certificates.jsonl').open('w') as out:
  for i,(p,r) in enumerate(zip(candidates,cpp)):
   assert int(r['candidate'])==p and r['d']==d and r['N']==N
   assert p>2*N and p%(2*d) not in (1,2*d-1)
   r['characteristic_gate']=True;r['congruence_gate']=True
   if r['status']=='composite_factor':
    factor=int(r['factor']);assert 1<factor<p and p%factor==0;r['accepted_status']='EXACT_COMPOSITE_FACTOR'
   else:
    assert r['status'] in ('nonzero_determinant','zero_determinant')
    t=time.monotonic();det=int(fmpz_mod_mat(rows,fmpz_mod_ctx(p)).det());r['flint_seconds']=time.monotonic()-t
    r['flint_determinant_mod_n']=str(det);assert det==int(r['determinant_mod_n'])
    r['independent_determinants_agree']=True
    r['accepted_status']='IF_PRIME_THEN_CONNECTED' if det else 'UNRESOLVED_ZERO_SELECTED_MINOR'
    if det:assert r['unit_pivots']==N-2 and math.gcd(det,p)==1
   r['source_sha256']=source_audit['source_sha256'];r['source_audit_file']=str(BASE/f'd{d}_source_audit.json')
   results.append(r);out.write(json.dumps(r)+'\n');out.flush()
   if (i+1)%10==0 or i+1==len(candidates):
    print(json.dumps({'stage':'cross_verified','d':d,'completed':i+1,'total':len(candidates),'elapsed':time.monotonic()-started}),flush=True)
    save(BASE/f'd{d}_verification_state.json',{'status':'CROSS_VERIFICATION_RUNNING','d':d,'completed':i+1,'total':len(candidates),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
 final={'status':'COMPLETE','d':d,'N':N,'source_sha256':source_audit['source_sha256'],'gmp_code_sha256':sha(BASE/'independent_modular_det.cpp'),'python_code_sha256':sha(Path(__file__)),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidate_count':len(candidates),'conditional_connected':sum(r['accepted_status']=='IF_PRIME_THEN_CONNECTED' for r in results),'composite':sum(r['accepted_status']=='EXACT_COMPOSITE_FACTOR' for r in results),'unresolved':sum(r['accepted_status']=='UNRESOLVED_ZERO_SELECTED_MINOR' for r in results),'independence':'GMP hand-written row-pivot Gaussian elimination versus FLINT determinant; distinct implementations, same audited integer source.','elapsed_seconds':time.monotonic()-started,'records':results}
 save(BASE/f'd{d}_complete_verification.json',final);save(BASE/f'd{d}_verification_state.json',{k:v for k,v in final.items() if k!='records'})
 print(json.dumps({k:v for k,v in final.items() if k!='records'},indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('d',type=int);ap.add_argument('--manifest',type=Path,required=True);args=ap.parse_args();verify_group(args.d,args.manifest.resolve())
