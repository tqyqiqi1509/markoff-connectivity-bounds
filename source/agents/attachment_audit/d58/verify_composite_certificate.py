#!/usr/bin/env python3
"""Independent source normalization and exact bounded-modular minor verification.
Generic composite d only; no codimension-one assertion is used.
"""
import argparse,datetime,functools,hashlib,json,math,sys,time
from pathlib import Path
from flint import fmpz_mat,nmod_mat
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2)+'\n');t.replace(p)
def primes(n):return [x for x in range(2,n+1) if all(x%q for q in range(2,math.isqrt(x)+1))]
def strip(x,ps):
 x=abs(x)
 if not x:return 0
 for q in ps:
  while x%q==0:x//=q
 return x
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out
def monic_div(a,b):
 assert b[-1]==1;out=[0]*(len(a)-len(b)+1);w=a[:]
 for i in reversed(range(len(out))):
  out[i]=w[i+len(b)-1]
  for j,v in enumerate(b):w[i+j]-=out[i]*v
 assert not any(w)
 return out
def h(r):
 out=[0]*(r//2+1)
 for j in range((r-1)//2+1):out[(r-1-2*j+int(r%2==0))//2]+=(-1)**j*math.comb(r-1-j,j)
 while len(out)>1 and out[-1]==0:out.pop()
 return out
def filt(d,n):
 qs=[q for q in primes(d) if d%q==0];assert len(qs)==2 and math.prod(qs)==d;rs=[]
 for q in qs:
  v=n
  while v%q==0:v//=q
  rs.append(v)
 return monic_div(mul(h(rs[0]),h(rs[1])),h(math.gcd(*rs)))
def had(rows):
 H=1
 for row in rows:
  v=sum(x*x for x in row);assert v>0;a=math.isqrt(v);H*=a+(a*a<v)
 return H

def main(d):
 st=time.monotonic();cp=HERE/f'composite_d{d}_certificate.json';c=json.load(cp.open());mp=HERE/c['matrix_file'];s=json.load(mp.open());N=c['nd'];bound=2*N;ps=primes(bound)
 assert c['status']=='certified_no_exception_p_above_bound' and c['rough_gcd']=='1' and c['certificate_method']=='maximal_minors'
 assert c['d']==s['d']==d and s['nd']==N and d == 46 and N==6*d and N>=3*d and s['kappa']==0 and c['matrix_sha256']==sha(mp)
 assert s['p_condition']==f'p>{bound} and p not congruent +/-1 mod {2*d}'
 params=[];expected=set()
 for n in range(d,N+1,d):
  g=filt(d,n);count=n//2-len(g)+1;params.append({'n':n,'m_count':count,'g_T':g});expected.update((n,m) for m in range(count))
 assert s['parameters']==params
 ap=HERE/s['filter_audit_file'];aud=json.load(ap.open());assert aud['d']==d and aud['nd']==N
 assert [(r['n'],r['m_count'],r['g_T']) for r in aud['records']]==[(p['n'],p['m_count'],p['g_T']) for p in params]
 for r in aud['records']:
  assert r['exact_filter_zero_set_verified'] and len(r['selected_beta_classes'])==r['m_count']
  classes=r['selected_beta_classes']+r['excluded_beta_classes'];assert sorted(x['representative_exponent'] for x in classes)==list(range(1,r['n']//2+1))
  for key,want in [('selected_beta_classes',True),('excluded_beta_classes',False)]:
   for x in r[key]:
    o=2*r['n']//math.gcd(x['representative_exponent'],2*r['n']);assert o==x['ordinary_order'] and math.lcm(o,2)==x['even_order'] and (x['even_order']%(2*d)==0)==want
 raw=HERE/s['raw_rows_file'];assert sha(raw)==s['raw_rows_sha256'];raw_by={};bins=[math.comb(2*i,i) for i in range(N+1)]
 for line in raw.open():
  assert line.endswith('\n');r=json.loads(line);key=r['n'],r['m'];assert r['d']==d and r['nd']==N and key in expected and key not in raw_by
  co=list(map(int,r['coefficients']));assert len(co)==key[0]+1 and co[0]==0 and co[1]==12*co[2] and sum(a*b for a,b in zip(co,bins))==0
  content=functools.reduce(math.gcd,co,0);rt=strip(content,ps);div=abs(content)//rt if content else 1;assert strip(div,ps)==1
  raw_by[key]=([co[j]//div if j<len(co) else 0 for j in range(3,N+1)],{'n':key[0],'m':key[1],'divided_small_content':str(div),'content_rough':str(rt)})
 assert set(raw_by)==expected
 ordered=[raw_by[k] for k in sorted(expected,reverse=True)];assert [x[0] for x in ordered]==s['rows'] and [x[1] for x in ordered]==s['row_metadata']
 assert c['shape']==[len(s['rows']),N-2]
 used=[];rough_gcd=0
 for idx,entry in enumerate(c['minor_certificates']):
  D=int(entry['det']);assert str(strip(D,ps))==entry['rough_det'];indices=entry['row_indices'];assert len(indices)==N-2 and len(set(indices))==N-2 and min(indices)>=0 and max(indices)<len(s['rows'])
  if not D:continue
  rows=[s['rows'][i] for i in indices];H=had(rows);assert abs(D)<=H;used.append({'index':idx,'D':D,'H':H,'matrix':fmpz_mat(rows),'checks':0});rough_gcd=math.gcd(rough_gcd,strip(D,ps))
  if rough_gcd==1:break
 assert rough_gcd==1 and used
 hp=HERE/f'composite_d{d}_modular_header.json';jp=HERE/f'composite_d{d}_modular_journal.jsonl';identity={'d':d,'N':N,'matrix_sha256':sha(mp),'raw_sha256':sha(raw),'certificate_sha256':sha(cp),'filter_audit_sha256':sha(ap),'code_sha256':sha(Path(__file__)),'minor_indices':[x['index'] for x in used],'Hadamard_bits':[x['H'].bit_length() for x in used]}
 if hp.exists():assert json.load(hp.open())==identity
 else:save(hp,identity)
 M=1;k=(1<<29)+1
 while k%6!=1:k+=2
 k-=6;count=0
 if jp.exists():
  for line in jp.open():
   assert line.endswith('\n');r=json.loads(line);q=r['q'];assert r['k']>k and q==r['k']*(1<<32)+1 and q<(1<<64) and r['k']%2 and pow(3,(q-1)//2,q)==q-1
   for u,v in zip(used,r['det_residues']):
    if M<=2*u['H']:assert v==u['D']%q;u['checks']+=1
    else:assert v is None
   M*=q;k=r['k'];count+=1
 with jp.open('a') as log:
  while M<=2*max(x['H'] for x in used):
   while True:
    k+=6;q=k*(1<<32)+1;assert q<(1<<64)
    if pow(3,(q-1)//2,q)==q-1:break
   vals=[]
   for u in used:
    if M<=2*u['H']:
     v=int(nmod_mat(u['matrix'],q).det());assert v==u['D']%q;vals.append(v);u['checks']+=1
    else:vals.append(None)
   log.write(json.dumps({'q':q,'k':k,'exponent':32,'witness':3,'det_residues':vals})+'\n');log.flush();M*=q;count+=1
   if count%100==0:print(json.dumps({'stage':'COMPOSITE_EXACT_MODULAR_MINOR_VERIFICATION','d':d,'q_count':count,'product_bits':M.bit_length(),'target_bits':(2*max(x['H'] for x in used)).bit_length(),'seconds':time.monotonic()-st}),flush=True)
 out={'status':'GENERAL_COMPOSITE_UNIFORM_CERTIFICATE_INDEPENDENTLY_VERIFIED',**identity,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'strict_p_bound':bound,'excluded_residues_modulus':2*d,'excluded_residues':[1,2*d-1],'raw_rows_and_normalization_independently_checked':True,'independent_filter_binomial_lcm_and_full_class_order_checks':True,'rough_gcd_of_verified_integer_minors':'1','modular_prime_count':count,'modular_journal_sha256':sha(jp),'modulus_product':str(M),'used_minor_proofs':[{'index':u['index'],'Hadamard_bound':str(u['H']),'checks':u['checks']} for u in used],'prime_proof':'q=k*2^32+1<2^64, k odd, 3^((q-1)/2)=-1 mod q forces every prime divisor to exceed sqrt(q).','determinant_proof':'Both integers lie within +/-H; exact congruences modulo distinct proved primes with product >2H force equality.','shared_arithmetic_backend':'FLINT; integer and modular determinant entry points differ. Independent source/certificate logic, not independent implementation of every base arithmetic operation.','scope':'For every prime p>2N outside +/-1 mod 2d, rank criterion applies; relies on separately audited general-d bootstrap and pure-Vieta bridge.','elapsed_seconds':time.monotonic()-st}
 save(HERE/f'composite_d{d}_independent_verification.json',out);print(json.dumps({k:v for k,v in out.items() if k not in ['modulus_product','used_minor_proofs']},indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('d',type=int);main(ap.parse_args().d)
