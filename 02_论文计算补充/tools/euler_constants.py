"""Directed Euler-factor checks for the constants used in the article."""
from pathlib import Path
from fractions import Fraction as F
import json,math,hashlib
from flint import arb,ctx
ctx.prec=256
def aa(x):
 x=F(x);return arb(x.numerator)/x.denominator
def primes(n):
 s=bytearray([1])*n;s[:2]=bytes(2)
 for q in range(2,math.isqrt(n-1)+1):
  if s[q]:s[q*q:n:q]=bytes(len(range(q*q,n,q)))
 return [q for q in range(2,n) if s[q]]
def pack(x):return {'upper':str(x.upper().fmpq()),'display':str(x)}
ps=primes(449);assert len(ps)==86 and 449**5>2**44
tau=arb(1);tau_rows=[]
for q in ps:
 step=(-aa(F(5,44))*arb(q).log()).exp();v=arb(1);best=arb(1);a=0
 while True:
  a+=1;v*=step;best=best.max((a+1)*v)
  if aa(F(a+2,a+1))*step<1:break
  assert a<1000
 tau*=best.upper();tau_rows.append({'q':q,'last_a':a,'local':pack(best)})
assert tau<2541933
delta=F(829,9180);eta=F(1,18);prod=arb(1);low=arb(1);B_rows=[]
assert (-aa(delta)*arb(401).log()).exp()*(1+(-aa(eta)*arb(401).log()).exp())<1
for q in primes(401):
 h=arb(q).log();r=(-aa(eta)*h).exp();step=(-aa(delta)*h).exp();c=1-arb(1)/q;inf=1+c*r/(1-r);B=arb(1);cur=arb(1);wp=arb(1);best=arb(1);a=0
 while True:
  a+=1;cur*=r;wp*=step;B+=c*cur;best=best.max(wp*B)
  if wp*inf<1:break
  assert a<1000
 prod*=best.upper()
 if q>=31:low*=best.upper()
 B_rows.append({'q':q,'tail_from':a,'local':pack(best)})
assert prod<159842 and low<662
assert 73**305<72**306
ratio=arb(1);ratio_rows=[]
for q in [q for q in primes(73) if q>=31]:
 E=1
 while q**E*(q-1)**306<q**306:E+=1
 h=arb(q).log();r=(-h/17).exp();s=(-h/18).exp();A=B=arb(1);rp=sp=arb(1);best=arb(1)
 for a in range(1,E+1):
  rp*=r;sp*=s;A+=rp;B+=(1-arb(1)/q)*sp;best=best.max(A/B)
 ratio*=best.upper();ratio_rows.append({'q':q,'tail_from':E,'local':pack(best)})
assert ratio<aa(F(26,25))
out={'status':'PASS_EULER_CONSTANTS_ALL_PRIME_AND_EXPONENT_TAILS','tau':pack(tau),'B_all':pack(prod),'B_support_ge31':pack(low),'low_count_ratio':pack(ratio),'tau_rows':tau_rows,'B_rows':B_rows,'ratio_rows':ratio_rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(Path(__file__).resolve().parent.parent/'certificates/EULER_CONSTANTS.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if not k.endswith('_rows')},indent=2))
