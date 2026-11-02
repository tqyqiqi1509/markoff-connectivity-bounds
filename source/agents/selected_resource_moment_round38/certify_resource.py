#!/usr/bin/env python3
"""Fresh finite resource proof; no prior producer/table imports.

Variance uses pairwise weighted differences; tilt uses an algebraic18th root.
All prime and exponent tails are independent of selected-variable enumeration.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,math,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb,ctx
ctx.prec=256
def aa(q):
 q=F(q);return arb(q.numerator)/q.denominator
def pp(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def prime(q):return q>=2 and all(q%d for d in range(2,math.isqrt(q)+1))
def moments(q,a):
 h=arb(q).log();r=1/arb(q).root(18);c=arb(q-1)/q
 ws=[arb(1)]+[c*r**j for j in range(1,a+1)];B=sum(ws,arb(0))
 mu=sum((j*w for j,w in enumerate(ws)),arb(0))/B
 var=sum((ws[i]*ws[j]*(j-i)**2 for i in range(a+1) for j in range(i+1,a+1)),arb(0))/B**2*h**2
 third=sum((w*abs(arb(j)-mu)**3 for j,w in enumerate(ws)),arb(0))/B*h**3
 return var,third
def main():
 start=time.monotonic();primes=[q for q in range(2,409) if prime(q)];assert len(primes)==79 and primes[-1]==401
 loglower=16*sum((F(1,3)**(2*j+1)/F(2*j+1) for j in range(8)),F(0))+2*sum((F(153,665)**(2*j+1)/F(2*j+1) for j in range(8)),F(0))
 assert loglower>6 and arb(409).log()>6
 rows=[];prime_rows=[];total=F(0)
 for q in primes:
  cap=F(0);winner=0;a=1
  while True:
   v,b=moments(q,a);f=3*v-b
   row={'q':q,'a':a,'v':pp(v),'b':pp(b),'resource':pp(f),'is_stop':bool(v>9)};rows.append(row)
   if v>9:
    assert f<0;stop=row;break
   upper=F(str(f.upper().fmpq()))
   if upper>cap:cap=upper;winner=a
   a+=1;assert a<=100,'Explicit unexpected-failure stop, never a mathematical exponent cutoff'
  total+=cap;prime_rows.append({'q':q,'stop_exponent':a,'stop_v_lower':str(v.lower().fmpq()),'cap_upper':str(cap),'maximizing_exponent_for_upper':winner})
 assert total<F(164)
 lam=F(43,1000);nu=F(61,10000);D=F(1693,200)
 floor=(D-nu*164)/(lam-3*nu);assert floor==F(5742,19) and floor>302
 out={'status':'RESOURCE_CONSTANT_CERTIFIED_PENDING_ROOT_ACCEPTANCE','scope':'All actual tilted prime-power variables, any finite subset of distinct primes, any actual nonnegative exponent. No selection condition or prior giant moment table required.','t':'3','C_upper':str(total),'C_display_only':str(aa(total)),'clean_C':'164','prime_tail_start':409,'log409_rational_lower':str(loglower),'log409_arb':pp(arb(409).log()),'finite_prime_count':len(primes),'finite_rows_including_stops':len(rows),'largest_stop':max(x['stop_exponent'] for x in prime_rows),'positive_prime_caps':sum(F(x['cap_upper'])>0 for x in prime_rows),'prime_rows':prime_rows,'moment_rows':rows,'lambda':'43/1000','nu':'61/10000','D':'1693/200','V_lower_from_clean_C':str(floor),'safe_V_floor':302,'strict_gap_at302':str(D-nu*164-(lam-3*nu)*302),'method':'Arb256 algebraic18th roots, weighted pairwise variance, direct absolute central third moment, exhaustive trial-division primes, fresh per-prime variance>=9 stops.','arithmetic_backend':'python-flint/Arb; same backend family as earlier proofs, no algorithmic import.','producer_sha256':sha(__file__),'seconds':time.monotonic()-start}
 (HERE/'RESOURCE_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['prime_rows','moment_rows']},indent=2),flush=True)
if __name__=='__main__':main()
