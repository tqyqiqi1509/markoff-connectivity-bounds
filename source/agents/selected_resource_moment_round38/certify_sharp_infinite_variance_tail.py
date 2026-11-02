#!/usr/bin/env python3
"""Same1600 split, stronger monotone BE infinite tail; no parameter sweep."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb,ctx
ctx.prec=256
def aa(x):x=F(x);return arb(x.numerator)/x.denominator
def pp(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def main():
 cert=json.loads((HERE/'GENERIC_RHO_TRACE_CERTIFICATE.json').read_text());rows=cert['rows']
 assert [r['V_interval'] for r in rows]==[[str(v),str(v+1)] for v in range(302,1600)]
 gamma=F(177,500);maximum=max(F(r['total_upper']['upper']) for r in rows);assert maximum<gamma
 v=1600;a=arb(v).sqrt()/18;rt2=arb(2).sqrt();rtpi=(2*arb.pi()).sqrt()
 def phi(w):w=aa(w);return (-w*w/2).exp()/rtpi
 def Phi(w):return (-(aa(w)/rt2)).erfc()/2
 def derivative(w):return phi(w)-a*Phi(w)
 lo=F(-10);hi=F(10);assert derivative(lo)>0 and derivative(hi)<0
 for _ in range(64):
  mid=(lo+hi)/2;val=derivative(mid)
  if val>0:lo=mid
  else:assert val<0;hi=mid
 w=arb(aa((lo+hi)/2),aa((hi-lo)/2))
 h=((-a*a/2-a*w).exp()*(-(w/rt2)).erfc()/2).upper()
 lam=F(43,1000);nu=F(61,10000);D=F(1693,200)
 assert 1600>3*D/lam
 error=aa(F(28,25)/nu*(lam*v-D))/arb(v).sqrt()**3
 tail=(h+error).upper();assert tail<aa(gamma)
 upper=max(maximum,F(str(tail.fmpq())))
 out={'status':'SHARPER_TRACE_COROLLARY_CERTIFIED_PENDING_ROOT_ACCEPTANCE','gamma':str(gamma),'finite_variance_domain':['302','1600'],'finite_boxes':len(rows),'finite_maximum':str(maximum),'finite_maximum_display':str(aa(maximum)),'tail_from':1600,'BE_critical_variance':str(3*D/lam),'normal_root_interval':[str(lo),str(hi)],'normal_H_upper':pp(h),'BE_error_upper':pp(error.upper()),'infinite_variance_tail_upper':pp(tail),'global_upper':str(upper),'global_display_only':str(aa(upper)),'argument':'H(sqrt(V)/18) is nonincreasing by Gaussian convolution; (lambda V-D)/V^(3/2) decreases for V>=3D/lambda. Thus fullBE(V)<=fullBE(1600) for allV>=1600. Keep same finite grid and genericrho.','bindings':{'GENERIC_RHO_TRACE_CERTIFICATE.json':sha(HERE/'GENERIC_RHO_TRACE_CERTIFICATE.json'),'RESOURCE_CERTIFICATE.json':sha(HERE/'RESOURCE_CERTIFICATE.json')},'producer_sha256':sha(__file__)}
 (HERE/'SHARP_INFINITE_VARIANCE_TAIL_COROLLARY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
