#!/usr/bin/env python3
"""Independent exact-integer interval replay. Standard library only; no Arb.
Every arithmetic operation rounds outwards to a fixed160-bit dyadic grid.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
HERE=Path(__file__).resolve().parent;S=1<<160
def ceildiv(n,d):return -((-n)//d)
class I:
 def __init__(self,lo,hi=None):self.lo=int(lo);self.hi=int(lo if hi is None else hi);assert self.lo<=self.hi
 @staticmethod
 def val(x):
  x=F(x);return I((x.numerator*S)//x.denominator,ceildiv(x.numerator*S,x.denominator))
 def __add__(self,o):
  o=o if isinstance(o,I) else I.val(o);return I(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+-o if isinstance(o,I) else self+(-F(o))
 def __mul__(self,o):
  o=o if isinstance(o,I) else I.val(o);v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
  return I(min(v)//S,ceildiv(max(v),S))
 __rmul__=__mul__
 def __truediv__(self,o):
  o=o if isinstance(o,I) else I.val(o);assert o.lo>0
  v=[F(a*S,b) for a in [self.lo,self.hi] for b in [o.lo,o.hi]];lo=min(v);hi=max(v)
  return I(lo.numerator//lo.denominator,ceildiv(hi.numerator,hi.denominator))
 def __pow__(self,k):
  assert isinstance(k,int) and k>=0;x=I.val(1);y=self
  while k:
   if k&1:x=x*y
   y=y*y;k//=2
  return x
 def abs(self):
  if self.lo>=0:return self
  if self.hi<=0:return -self
  return I(0,max(-self.lo,self.hi))
 def pack(self):return {'lower':str(F(self.lo,S)),'upper':str(F(self.hi,S))}
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def atanh_log(z):
 assert 0<=z<=F(1,3);n=48
 low=2*sum((z**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
 hi=low+2*z**(2*n+1)/(F(2*n+1)*(1-z*z))
 return low,hi
L2=atanh_log(F(1,3))
def logq(q):
 k=q.bit_length()-1;two=1<<k;l,u=atanh_log(F(q-two,q+two))
 lo=k*L2[0]+l;hi=k*L2[1]+u
 return I((lo.numerator*S)//lo.denominator,ceildiv(hi.numerator*S,hi.denominator))
def tilt(q):
 lo=0;hi=S
 while hi-lo>1:
  m=(lo+hi)//2
  if q*m**18<=S**18:lo=m
  else:hi=m
 assert q*lo**18<=S**18<q*(lo+1)**18
 return I(lo,lo+1)
def moments(q,a):
 r=tilt(q);h=logq(q);c=I.val(F(q-1,q));weights=[I.val(1)]+[c*r**j for j in range(1,a+1)]
 B=sum(weights,I.val(0));mu=sum((w*j for j,w in enumerate(weights)),I.val(0))/B
 second=sum((w*(j*j) for j,w in enumerate(weights)),I.val(0))/B
 v=(second-mu**2)*h**2
 b=sum((w*(I.val(j)-mu).abs()**3 for j,w in enumerate(weights)),I.val(0))/B*h**3
 assert v.lo>0 and b.lo>0
 return v,b
def main():
 st=time.monotonic();p=HERE/'RESOURCE_CERTIFICATE.json';d=json.loads(p.read_text());assert d['t']=='3' and d['clean_C']=='164'
 sieve=[True]*409;sieve[0]=sieve[1]=False
 for k in range(2,21):
  if sieve[k]:
   for j in range(k*k,409,k):sieve[j]=False
 primes=[q for q in range(409) if sieve[q]];assert primes==[r['q'] for r in d['prime_rows']] and len(primes)==79
 assert logq(409).lo>6*S
 oldrows={(r['q'],r['a']):r for r in d['moment_rows']};assert len(oldrows)==194
 newrows=[];caps=[];total=0;seen=set()
 for q,pr in zip(primes,d['prime_rows']):
  best=0;winner=0;stop=pr['stop_exponent']
  for a in range(1,stop+1):
   v,b=moments(q,a);f=3*v-b;seen.add((q,a));old=oldrows[q,a]
   for key,val in [('v',v),('b',b),('resource',f)]:
    assert F(old[key]['lower'])<=F(val.hi,S) and F(val.lo,S)<=F(old[key]['upper']),(q,a,key,'disjoint enclosures')
   if a==stop:assert v.lo>9*S and f.hi<0
   else:
    assert v.hi<9*S
    if f.hi>best:best=f.hi;winner=a
   newrows.append({'q':q,'a':a,'v':v.pack(),'b':b.pack(),'resource':f.pack(),'is_stop':a==stop})
  total+=best;caps.append({'q':q,'cap_upper':str(F(best,S)),'maximizing_exponent_for_upper':winner,'stop_exponent':stop})
 assert seen==set(oldrows) and total<164*S
 assert F(d['C_upper'])<164 and F(1693,200)-F(61,10000)*164-F(302)*(F(43,1000)-3*F(61,10000))==F(13,2500)
 out={'status':'PASS_INDEPENDENT_PURE_RATIONAL_RESOURCE_REPLAY','arithmetic_backend':'Python exact integers/Fractions only; no Arb, FLINT, floating-point gates or producer imports.','precision_bits':160,'log_method':'48-term positive atanh series with geometric remainder after power-of-two reduction','moment_method':'Raw first/second moments and absolute centered third moment; differs from producer pairwise variance.','prime_count':79,'moment_rows':194,'C_upper':str(F(total,S)),'C_display_only':float(F(total,S)),'clean_C':164,'V_lower_clean':'5742/19','safe_V_floor':302,'prime_caps':caps,'rows':newrows,'source_sha256':sha(p),'verifier_sha256':sha(__file__),'seconds':time.monotonic()-st,'cross_check_scope':'Fresh rigorous C<164 and all prime/exponent tails. Original higher-precision Arb enclosures overlap the independent dyadic enclosures; no claim to prove every original tighter endpoint by containment.'}
 (HERE/'RATIONAL_RESOURCE_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['prime_caps','rows']},indent=2),flush=True)
if __name__=='__main__':main()
