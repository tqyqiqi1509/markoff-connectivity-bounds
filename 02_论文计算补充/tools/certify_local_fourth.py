#!/usr/bin/env python3
"""All-q/all-a local fourth moment lemma. No numerical endpoint claim.

For selected actual divisor factors (43/1000)v-(61/10000)b>0,
certifies E[X^4] <= 8 E|X|^3.  All a=1 and q>=2^15 use
analytic lemmas recorded in PROOF_EDGEWORTH.md.  The finite part is fresh.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,sys,platform
from flint import arb,ctx
ctx.prec=192
HERE=Path(__file__).resolve().parent
LAM=F(43,1000);NU=F(61,10000);C=LAM/NU;RHO=C*C;Q=2**15

def aa(q):
 q=F(q);return arb(q.numerator)/q.denominator

def pp(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq())}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def primes(limit):
 s=bytearray(b'\1')*limit;s[0:2]=b'\0\0'
 for p in range(2,int(limit**.5)+1):
  if s[p]:s[p*p:limit:p]=b'\0'*len(range(p*p,limit,p))
 return [p for p in range(2,limit) if s[p]]

def main():
 st=time.monotonic();ps=primes(Q);rows=[];stops=[];selected=excluded=0;max_ratio=F(0);worst=None
 # Entire q>=Q domain: lattice support gives h<2C if selected.
 # Then r>9/20 by a rational two-term atanh lower bound for log(20/9).
 z=F(11,29);assert 2*(z+z**3/3)>C/9
 assert F(3,5)**6>F(1,32)
 # c_q>=999/1000; polynomial numerator lower bound >=49/100 denominator.
 r=F(9,20);cq=F(999,1000)
 poly=cq*r+4*cq*r*r+cq*cq*r**3-F(49,100)*(1+r+r*r)**2
 derivative_lower=F(19,1000)+F(5052,1000)*F(9,20)-F(196,100)*F(3,5)**3
 assert poly>0 and derivative_lower>0
 log2lower=2*sum((F(1,3)**(2*j+1)/F(2*j+1) for j in range(3)),F(0))
 assert log2lower>F(69,100)
 assert F(49,100)*(15*F(69,100))**2>RHO
 # Every two-point distribution: mu4*v/b^2 <=9/8, then mu4/b <9C/8<8.
 assert 9*C/8<8
 for q in ps:
  h=arb(q).log();r=(-h/18).exp();cq=1-arb(1)/q
  ws=[arb(1)];B=arb(1);s1=arb(0);s2=arb(0);cur=arb(1)
  for a in range(1,257):
   cur*=r;w=cq*cur;ws.append(w);B+=w;s1+=a*w;s2+=a*a*w
   mu=s1/B;v=(s2/B-mu*mu)*h*h
   if v>aa(RHO):
    stops.append({'q':q,'a':a,'variance':pp(v)});break
   if a==1:continue
   b=sum((wj*abs(arb(j)-mu)**3 for j,wj in enumerate(ws)),arb(0))*h**3/B
   u4=sum((wj*(arb(j)-mu)**4 for j,wj in enumerate(ws)),arb(0))*h**4/B
   sel=aa(LAM)*v-aa(NU)*b
   if sel<0:
    status='excluded_by_selection';excluded+=1
   else:
    assert u4<8*b, (q,a,str(u4-8*b),str(sel))
    status='selected_or_uncertain_fourth_verified';selected+=1
    ratio=F(str((u4/b).upper().fmpq()))
    if ratio>max_ratio:max_ratio=ratio;worst={'q':q,'a':a,'ratio':str(ratio),'display_only':str(u4/b)}
   rows.append({'q':q,'a':a,'variance':pp(v),'abs_third':pp(b),'fourth':pp(u4),'selection':pp(sel),'status':status})
  else:raise AssertionError(('missing_variance_stop',q))
 assert len(stops)==len(ps)
 out={'status':'PASS_LOCAL_FOURTH_ALL_PRIMES_ALL_EXPONENTS_WITH_ANALYTIC_TAILS','statement':'selected actual local factors have E[X^4] < 8 E|X|^3','parameters':{'lambda':str(LAM),'nu':str(NU),'c':str(C),'rho':str(RHO),'q_cutoff':Q},'prime_count':len(ps),'finite_a_ge_2_rows':len(rows),'selected_or_uncertain_verified':selected,'excluded_rows':excluded,'maximum_stop_exponent':max(x['a'] for x in stops),'finite_maximum_fourth_over_abs_third':worst,'analytic_exact_checks':{'two_point_ratio_bound':'9/8','two_point_fourth_over_abs_third':str(9*C/8),'log20over9_lower':str(2*(z+z**3/3)),'a2_variance_rational_poly_at_9over20':str(poly),'a2_poly_derivative_lower':str(derivative_lower),'large_q_variance_lower':str(F(49,100)*(15*F(69,100))**2)},'rows':rows,'stops':stops,'source_sha256':sha(__file__),'prime_list_sha256':hashlib.sha256(','.join(map(str,ps)).encode()).hexdigest(),'environment':{'python':sys.version,'platform':platform.platform(),'arb_precision_bits':ctx.prec},'seconds':time.monotonic()-st}
 (HERE/'LOCAL_FOURTH_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k not in ['rows','stops']},indent=2),flush=True)
if __name__=='__main__':main()
