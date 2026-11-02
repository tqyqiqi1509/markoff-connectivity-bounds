#!/usr/bin/env python3
"""Combine disjoint raw prefix/tail sources and divide only <=8d-smooth content."""
import argparse,datetime,functools,hashlib,json,math,sys,time
from pathlib import Path
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent

def primes(n):
 a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
 for q in range(2,math.isqrt(n)+1):
  if a[q]:a[q*q:n+1:q]=b'\0'*((n-q*q)//q+1)
 return [i for i in range(2,n+1) if a[i]]

def strip(n,ps):
 n=abs(n)
 if not n:return 0
 for p in ps:
  while n%p==0:n//=p
 return n

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def build(d):
 N=4*d;assert d>=5 and all(d%q for q in range(2,math.isqrt(d)+1))
 ps=primes(2*N);expected={(k*d,m) for k in range(1,5) for m in range(k*(d-1)//2)}
 indexed={};raw_provenance=[];total_states=total_adds=0;row_seconds=0.;start=time.monotonic();binoms=[math.comb(2*i,i) for i in range(N+1)]
 for part in ['prefix','tail']:
  path=BASE/f'd{d}_{part}_rows.jsonl'
  for line in path.open():
   assert line.endswith('\n');r=json.loads(line);assert r['d']==d and r['nd']==N
   key=r['n'],r['m'];assert key in expected and key not in indexed
   assert (part=='prefix')==(key[0]<N or key[1]<=d)
   co=list(map(int,r['coefficients']));assert len(co)==r['n']+1
   assert co[0]==0 and co[1]==12*co[2] and sum(c*b for c,b in zip(co,binoms))==0
   content=functools.reduce(math.gcd,co,0);rough=strip(content,ps);div=abs(content)//rough if content else 1
   assert div>=1 and strip(div,ps)==1 and all(c%div==0 for c in co)
   row=[(co[i] if i<len(co) else 0)//div for i in range(3,N+1)]
   indexed[key]=(row,{'n':key[0],'m':key[1],'divided_small_content':str(div),'content_rough':str(rough),'raw_partition':part})
   total_states+=r['nonzero_states'];total_adds+=r['additions'];row_seconds+=r['seconds']
  raw_provenance.append({'file':path.name,'sha256':sha(path),'bytes':path.stat().st_size})
 assert set(indexed)==expected
 rows=[];meta=[]
 for key in sorted(expected,reverse=True):
  row,r=indexed[key];rows.append(row);meta.append(r)
 gs=[[1],[0,1],[-1,1],[0,-2,1]]
 obj={'d':d,'nd':N,'kappa':0,'p_condition':f'p>{2*N} and p not congruent +/-1 mod {2*d}','rows':rows,'row_metadata':meta,'parameters':[{'n':k*d,'m_count':k*(d-1)//2,'g_T':gs[k-1]} for k in range(1,5)],'raw_sources':raw_provenance,'engine':'Audited C++ GMP descending half-degree Phi propagation, unchanged mathematical reduction','generator_versions':{'prefix':sha(BASE/'phi_stream_selected.cpp'),'tail':sha(BASE/'phi_stream_partitioned.cpp')},'all_three_integer_moments_verified':True,'all_allowed_rows_present_once':True}
 out=BASE/f'full_d{d}_matrix.json'
 with out.open('w') as f:json.dump(obj,f,separators=(',',':'))
 rec={'status':'FULL_SOURCE_AND_NORMALIZATION_PASSED','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'d':d,'N':N,'rows':len(rows),'columns':N-2,'matrix_file':out.name,'matrix_sha256':sha(out),'matrix_bytes':out.stat().st_size,'raw_sources':raw_provenance,'nonzero_states_total':total_states,'integer_additions_total':total_adds,'sum_row_seconds':row_seconds,'build_seconds':time.monotonic()-start,'normalization_only_primes_le':2*N}
 (BASE/f'd{d}_full_matrix_build.json').write_text(json.dumps(rec,indent=2)+'\n')
 print(json.dumps(rec,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('d',type=int);build(ap.parse_args().d)
