#!/usr/bin/env python3
"""Independent certificate logic with exact Hadamard + proven-prime modular determinants.
The producer's fmpz integer determinant is not rerun. This verifier uses nmod
field determinants and an explicit integer bound. FLINT arithmetic is shared;
source reconstruction, certificate logic, and prime proofs are independent.
"""
import argparse,datetime,functools,hashlib,json,math,sys,time,fcntl,signal
import multiprocessing as mp
from fractions import Fraction
from pathlib import Path
from flint import fmpz_mat,nmod_mat,ctx
sys.set_int_max_str_digits(0)
SERIAL_VERIFIER_SHA='b7f0ec03ab61b29e2976172722b8af022d2f43c332eed92db7989cfb62c5c424'
IMPLEMENTATION_ID='proven_proth_parallel_nmod_v1'
_MATRICES=None
_DS=None

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while True:
   chunk=f.read(1<<20)
   if not chunk:break
   h.update(chunk)
 return h.hexdigest()

def small_primes(n):
 ps=[]
 for x in range(2,n+1):
  if all(x%q for q in ps if q*q<=x):ps.append(x)
 return ps

def rough(x,ps):
 x=abs(x)
 if x==0:return 0
 for q in ps:
  while x%q==0:x//=q
 return x

def write_json(path,obj):
 t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(path)

def hadamard(rows):
 h=1
 for row in rows:
  n=sum(a*a for a in row);assert n>0
  s=math.isqrt(n);r=s+(s*s<n);assert (r-1)**2<n<=r*r
  h*=r
 return h

def proth_certificate(q,k,a=3):
 return k>0 and k%2==1 and q==k*(1<<32)+1 and q<(1<<64) and pow(a,(q-1)//2,q)==q-1


def planned_tasks(product,lastk,Hs,offset=0):
 tasks=[];k=lastk;M=product
 while M<=2*max(Hs.values()):
  while True:
   k+=6;q=k*(1<<32)+1;assert q<(1<<64)
   if proth_certificate(q,k):break
  tasks.append((offset+len(tasks)+1,q,k,M<=2*Hs['lower'],M<=2*Hs['top']))
  M*=q
 return tasks,M

def initialize_worker():
 # Only the parent writes. Fork inherits completed, read-only matrices.
 signal.signal(signal.SIGINT,signal.SIG_IGN)
 signal.signal(signal.SIGTERM,signal.SIG_DFL)
 ctx.threads=1

def modular_worker(task):
 index,q,k,need_lower,need_top=task
 assert proth_certificate(q,k)
 rec={'q':q,'k':k,'proth_exponent':32,'witness':3}
 for name,required in [('lower',need_lower),('top',need_top)]:
  if required:
   value=int(nmod_mat(_MATRICES[name],q).det())
   assert value==_DS[name]%q,(name,q,value,_DS[name]%q)
   rec[name+'_det_mod_q']=value
  else:rec[name+'_det_mod_q']=None
 return index,rec

def interrupted(signum,frame):
 # Raising inside the Pool context terminates and joins all workers.
 raise SystemExit(128+signum)

def verify(d,source_dir,output_dir,workers):
 global _MATRICES,_DS
 BASE=source_dir.resolve();OUT=output_dir.resolve()
 assert BASE.is_dir() and OUT!=BASE and 1<=workers<=8
 assert digest(BASE/'verify_universal_modular.py')==SERIAL_VERIFIER_SHA
 OUT.mkdir(parents=True,exist_ok=True)
 lock=(OUT/'parallel_verifier.lock').open('a+')
 fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
 signal.signal(signal.SIGINT,interrupted);signal.signal(signal.SIGTERM,interrupted)
 code_sha=digest(Path(__file__))
 started=time.monotonic();cp=BASE/f'd{d}_universal_certificate.json';c=json.load(cp.open());src=BASE/c['source_file'];source_sha=digest(src)
 assert source_sha==c['source_sha256'] and c['d']==d and c['N']==4*d and c['bound']==8*d
 assert c['status']=='CERTIFIED_ALL_ELIGIBLE_PRIMES' and all(d%q for q in range(2,math.isqrt(d)+1))
 s=json.load(src.open());N=4*d;ps=small_primes(8*d)
 assert s['d']==d and s['nd']==N and s['kappa']==0
 expected=[(k*d,m) for k in range(1,5) for m in range(k*(d-1)//2)];assert len(expected)==5*d-5
 keys=[(m['n'],m['m']) for m in s['row_metadata']];assert keys==sorted(expected,reverse=True)
 assert len(s['rows'])==len(keys) and all(len(row)==N-2 for row in s['rows'])
 where={key:i for i,key in enumerate(keys)};seen=set();bins=[math.comb(2*i,i) for i in range(N+1)]
 for rsrc in s['raw_sources']:
  raw=BASE/rsrc['file'];assert digest(raw)==rsrc['sha256'] and raw.stat().st_size==rsrc['bytes']
  for line in raw.open():
   r=json.loads(line);key=r['n'],r['m'];assert r['d']==d and r['nd']==N and key in where and key not in seen;seen.add(key)
   co=[int(v) for v in r['coefficients']];assert len(co)==key[0]+1
   assert co[0]==0 and co[1]==12*co[2] and sum(a*b for a,b in zip(co,bins))==0
   i=where[key];meta=s['row_metadata'][i];content=functools.reduce(math.gcd,co,0);rt=rough(content,ps);div=abs(content)//rt if content else 1
   assert str(div)==meta['divided_small_content'] and str(rt)==meta['content_rough'] and rough(div,ps)==1
   assert all(v%div==0 for v in co)
   row=[(co[j]//div if j<len(co) else 0) for j in range(3,N+1)]
   assert row==s['rows'][i] and all(v==0 for v in row[key[0]-2:])
 assert seen==set(expected)
 low=[row[:3*d-2] for row,meta in zip(s['rows'],s['row_metadata']) if meta['n']<=3*d]
 lc=c['lower_certificate'];assert lc['source_sha256']==source_sha and lc['bound']==8*d
 w=[int(v) for v in lc['primitive_kernel']];L=int(lc['last_coordinate']);assert len(w)==3*d-2 and w[-1]==L>0 and functools.reduce(math.gcd,w,0)==1
 assert all(sum(a*b for a,b in zip(row,w))==0 for row in low)
 Dlo=int(lc['first_square_minor_det']);assert Dlo%L==0
 assert str(abs(Dlo)//L)==lc['gcd_all_maximal_minors'] and rough(abs(Dlo)//L,ps)==1 and lc['gcd_rough_above_bound']=='1'
 top=[row for row,meta in zip(s['rows'],s['row_metadata']) if meta['n']==N];r=len(top);assert r==2*d-2==c['top_block_size']
 Ahi=[row[-r:][::-1] for row in top];Dhi=int(c['top_block_determinant']);assert Dhi!=0 and rough(Dhi,ps)==1
 rhs=[sum(a*b for a,b in zip(row[:3*d-2],w)) for row in top];assert list(map(str,rhs))==c['top_rhs']
 z=[Fraction(v) for v in c['top_solution']];assert len(z)==r
 den=1
 for v in z:assert rough(v.denominator,ps)==1;den=math.lcm(den,v.denominator)
 zi=[v.numerator*(den//v.denominator) for v in z]
 assert all(sum(a*b for a,b in zip(row,zi))==den*b for row,b in zip(Ahi,rhs))
 assert c['residual_indices']==list(range(d,r)) and c['rough_gcd']=='1'
 assert functools.reduce(math.gcd,(rough(v.numerator,ps) for v in z[d:]),0)==1
 Alo=[row[:-1] for row in low];Hlo=hadamard(Alo);Hhi=hadamard(Ahi);assert abs(Dlo)<=Hlo and abs(Dhi)<=Hhi
 matrices={'lower':fmpz_mat(Alo),'top':fmpz_mat(Ahi)};Ds={'lower':Dlo,'top':Dhi};Hs={'lower':Hlo,'top':Hhi}
 journal=OUT/f'd{d}_modular_det_journal.jsonl';header=OUT/f'd{d}_modular_det_header.json'
 identity={'implementation_id':IMPLEMENTATION_ID,'serial_predecessor_sha256':SERIAL_VERIFIER_SHA,'source_sha256':source_sha,'certificate_sha256':digest(cp),'code_sha256':digest(Path(__file__)),'d':d,'lower_H_bits':Hlo.bit_length(),'top_H_bits':Hhi.bit_length()}
 if header.exists():assert json.load(header.open())==identity
 else:write_json(header,identity)
 M=1;counts={'lower':0,'top':0};bounds_done={};lastk=(1<<29)+1
 while lastk%6!=1:lastk+=2
 lastk-=6;reused=0
 if journal.exists():
  for line in journal.open():
   assert line.endswith('\n'), 'partial journal must be preserved and repaired explicitly before resume'
   item=json.loads(line);q=item['q'];k=item['k'];assert k>lastk and proth_certificate(q,k,item['witness'])
   for name in ['lower','top']:
    required=M<=2*Hs[name]
    if required:
     assert item[name+'_det_mod_q']==Ds[name]%q;counts[name]+=1
    else:assert item[name+'_det_mod_q'] is None
   M*=q;lastk=k;reused+=1
   for name in ['lower','top']:
    if name not in bounds_done and M>2*Hs[name]:bounds_done[name]=str(M)
 print(json.dumps({'stage':'ALL_SOURCE_AND_ALGEBRA_IDENTITIES_PASSED','d':d,'H_bits':{n:h.bit_length() for n,h in Hs.items()},'reused_exact_modular_records':reused}),flush=True)
 state=OUT/f'd{d}_independent_verification_state.json';tested=reused;k=lastk
 _MATRICES=matrices;_DS=Ds
 tasks,planned_product=planned_tasks(M,lastk,Hs,tested)
 print(json.dumps({'stage':'PARALLEL_TASK_PLAN_COMPLETE','d':d,'workers':workers,'fresh_tasks':len(tasks),'reused_records':reused,'planned_product_bits':planned_product.bit_length()}),flush=True)
 with journal.open('a') as out:
  if tasks:
   with mp.get_context('fork').Pool(workers,initializer=initialize_worker) as pool:
    for task,result in zip(tasks,pool.imap(modular_worker,tasks,chunksize=1)):
     index,rec=result;q=rec['q'];k=rec['k']
     assert index==tested+1==task[0] and (q,k)==task[1:3]
     assert proth_certificate(q,k,rec['witness']) and rec['proth_exponent']==32
     for name,planned_required in [('lower',task[3]),('top',task[4])]:
      required=M<=2*Hs[name];assert required==planned_required
      value=rec[name+'_det_mod_q']
      if required:
       assert type(value) is int and value==Ds[name]%q
       counts[name]+=1
      else:assert value is None
     out.write(json.dumps(rec)+'\n');out.flush();M*=q;tested+=1
     for name in ['lower','top']:
      if name not in bounds_done and M>2*Hs[name]:bounds_done[name]=str(M)
     if tested%100==0 or M>2*max(Hlo,Hhi):
      progress={'stage':'EXACT_PARALLEL_MODULAR_DETERMINANT_VERIFICATION','d':d,'workers':workers,'moduli_count':tested,'M_bits':M.bit_length(),'target_bits':(2*max(Hlo,Hhi)).bit_length(),'checks_per_matrix':counts,'elapsed_seconds':time.monotonic()-started,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};write_json(state,progress);print(json.dumps(progress),flush=True)
 assert M==planned_product and tested==reused+len(tasks)
 assert digest(src)==source_sha and digest(cp)==identity['certificate_sha256']
 assert digest(Path(__file__))==code_sha and digest(BASE/'verify_universal_modular.py')==SERIAL_VERIFIER_SHA
 for rsrc in s['raw_sources']:assert digest(BASE/rsrc['file'])==rsrc['sha256']
 assert set(bounds_done)=={'lower','top'}
 proof={'implementation_id':IMPLEMENTATION_ID,'parallel_workers':workers,'serial_predecessor_sha256':SERIAL_VERIFIER_SHA,'status':'UNIVERSAL_CERTIFICATE_INDEPENDENTLY_VERIFIED','d':d,'N':N,'bound':8*d,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':source_sha,'certificate_sha256':digest(cp),'verifier_sha256':digest(Path(__file__)),'all_raw_rows_reconstructed_and_normalizations_checked':True,'primitive_lower_kernel_verified':True,'lower_maximal_minor_gcd_smooth':True,'top_integer_system_verified_after_common_denominator':True,'residual_rough_gcd_one':True,'integer_determinants_verified_by_explicit_modular_bound':True,'shared_arithmetic_backend':'FLINT; producer fmpz_mat integer determinant, verifier nmod_mat field determinants. Separate certificate logic and integer bounds; not advertised as independent implementations of every arithmetic operation.','proth_prime_proof':'q=k*2^32+1<2^64, k odd, a^((q-1)/2)=-1 mod q implies every prime divisor r of q is 1 mod 2^32 and hence exceeds sqrt(q), so q is prime.','journal_file':journal.name,'journal_sha256':digest(journal),'proven_moduli_count':tested,'checks_per_matrix':counts,'Hadamard_bounds':{n:str(h) for n,h in Hs.items()},'first_sufficient_modulus_products':bounds_done,'elapsed_seconds':time.monotonic()-started}
 write_json(OUT/f'd{d}_independent_universal_verification.json',proof);write_json(state,{k:v for k,v in proof.items() if k not in ['Hadamard_bounds','first_sufficient_modulus_products']});print(json.dumps({k:v for k,v in proof.items() if k not in ['Hadamard_bounds','first_sufficient_modulus_products']},indent=2),flush=True)
 write_json(OUT/'PARALLEL_STATE.json',{'stage':'PARALLEL_COMPLETE_PENDING_ROOT_ACCEPTANCE','d':d,'implementation_id':IMPLEMENTATION_ID,'verifier_sha256':code_sha,'source_sha256':source_sha,'certificate_sha256':digest(cp),'verification_sha256':digest(OUT/f'd{d}_independent_universal_verification.json'),'journal_sha256':digest(journal),'moduli':tested,'workers':workers})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('d',type=int);ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--workers',type=int,default=4);a=ap.parse_args();verify(a.d,a.source_dir,a.output_dir,a.workers)
