#!/usr/bin/env python3
"""Exact all-eligible-prime rank certificate; no prior directory writes.
Reuses the audited codimension-one argument, but computes the top determinant
directly, without importing the symbolic top-half determinant formula.
"""
import argparse,datetime,functools,hashlib,json,math,sys,time
from pathlib import Path
from flint import fmpz_mat,fmpq_mat
from build_full_matrix import primes,strip,sha
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent

def save(path,obj):
 t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(path)

def run(d):
 src=BASE/f'full_d{d}_matrix.json';fingerprint=sha(src);source=json.load(src.open());N=4*d;bound=8*d;ps=primes(bound);st=time.monotonic()
 assert source['d']==d and source['nd']==N and len(source['rows'])==5*d-5
 statepath=BASE/f'd{d}_universal_state.json'
 def stage(name,**kw):
  obj={'d':d,'stage':name,'source_sha256':fingerprint,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-st,**kw};save(statepath,obj);print(json.dumps(obj),flush=True)
 lowerrows=[row[:3*d-2] for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']<=3*d]
 assert len(lowerrows)==3*d-3
 lp=BASE/f'd{d}_lower_kernel_checkpoint.json'
 low=json.load(lp.open()) if lp.exists() else {}
 if low.get('source_sha256')!=fingerprint:low={'d':d,'bound':bound,'source_sha256':fingerprint}
 if 'primitive_kernel' not in low:
  stage('LOWER_RATIONAL_SOLVE_STARTED');t=time.monotonic()
  A=fmpz_mat([row[:-1] for row in lowerrows]);b=fmpz_mat([[-row[-1]] for row in lowerrows]);x=A.solve(b)
  assert fmpq_mat(A)*x==fmpq_mat(b)
  L=1
  for i in range(x.nrows()):L=math.lcm(L,int(x[i,0].denominator))
  w=[int(x[i,0].numerator)*(L//int(x[i,0].denominator)) for i in range(x.nrows())]+[L]
  assert functools.reduce(math.gcd,w,0)==1 and all(sum(a*b for a,b in zip(row,w))==0 for row in lowerrows)
  low.update({'primitive_kernel':list(map(str,w)),'last_coordinate':str(L),'solve_seconds':time.monotonic()-t});save(lp,low)
 else:
  w=list(map(int,low['primitive_kernel']));L=int(low['last_coordinate']);assert w[-1]==L>0 and functools.reduce(math.gcd,w,0)==1
  assert all(sum(a*b for a,b in zip(row,w))==0 for row in lowerrows)
 if 'first_square_minor_det' not in low:
  stage('LOWER_INTEGER_DETERMINANT_STARTED');t=time.monotonic();D=int(fmpz_mat([row[:-1] for row in lowerrows]).det())
  assert D%L==0;g=abs(D)//L
  low.update({'first_square_minor_det':str(D),'gcd_all_maximal_minors':str(g),'gcd_rough_above_bound':str(strip(g,ps)),'determinant_seconds':time.monotonic()-t});save(lp,low)
 if low['gcd_rough_above_bound']!='1':
  stage('STOP_LOWER_BLOCK_OBSTRUCTION',rough_gcd=low['gcd_rough_above_bound']);return
 stage('LOWER_BLOCK_CERTIFIED',lower_det_digits=len(low['first_square_minor_det']),last_coordinate_digits=len(low['last_coordinate']))
 top=[row for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']==N];r=len(top);assert r==2*d-2
 A=fmpz_mat([row[-r:][::-1] for row in top]);rhs=[sum(a*b for a,b in zip(row[:3*d-2],w)) for row in top]
 tp=BASE/f'd{d}_top_checkpoint.json';up=json.load(tp.open()) if tp.exists() else {}
 if up.get('source_sha256')!=fingerprint:up={'source_sha256':fingerprint,'d':d,'r':r}
 if 'determinant' not in up:
  stage('TOP_INTEGER_DETERMINANT_STARTED');t=time.monotonic();Dtop=int(A.det());assert Dtop!=0 and strip(Dtop,ps)==1
  up.update({'determinant':str(Dtop),'determinant_seconds':time.monotonic()-t});save(tp,up)
 if 'solution' not in up:
  stage('TOP_RATIONAL_SOLVE_STARTED');t=time.monotonic();z=A.solve(fmpz_mat([[a] for a in rhs]));assert fmpq_mat(A)*z==fmpq_mat([[a] for a in rhs])
  for i in range(r):assert strip(int(z[i,0].denominator),ps)==1
  up.update({'rhs':list(map(str,rhs)),'solution':[str(z[i,0]) for i in range(r)],'solve_seconds':time.monotonic()-t});save(tp,up)
 from fractions import Fraction
 z=[Fraction(v) for v in up['solution']];g=functools.reduce(math.gcd,(strip(v.numerator,ps) for v in z[d:]),0)
 out={'d':d,'N':N,'bound':bound,'source_file':src.name,'source_sha256':fingerprint,'status':'CERTIFIED_ALL_ELIGIBLE_PRIMES' if g==1 else 'RESIDUAL_OBSTRUCTION','lower_certificate':low,'top_block_size':r,'top_block_determinant':up['determinant'],'top_rhs':up['rhs'],'top_solution':up['solution'],'residual_indices':list(range(d,r)),'rough_gcd':str(g),'top_determinant_method':'direct FLINT integer determinant; no top-half formula used','p_conditions':{'strict_lower_bound':bound,'excluded_residues_modulus':2*d,'excluded_residues':[1,2*d-1]},'conclusion':'For every prime p>8d with p not congruent to +/-1 mod 2d, the nonorigin normalized Markoff surface pure-Vieta graph is connected, conditional on the audited rank-to-graph theorem.','producer_code_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-st,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 save(BASE/f'd{d}_universal_certificate.json',out);stage(out['status'],rough_gcd=str(g))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('d',type=int);run(ap.parse_args().d)
