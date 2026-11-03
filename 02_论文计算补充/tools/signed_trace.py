#!/usr/bin/env python3
"""Rigorous signed-first-Edgeworth / CL transport. Fixed old budget, Omega=32.

New fourth-moment bound and positive third cumulants retain the leading phase.
Reference is a signed Gaussian derivative measure, not a probability measure.
Every frequency rectangle and variance interval is enclosed by Arb balls.
"""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,sys,platform
from flint import arb,ctx
ctx.prec=192
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
TABLE=ROOT/'certificates/FREQUENCY_ENVELOPE.json'
TSHA=hashlib.sha256(TABLE.read_bytes()).hexdigest()
ETA=F(1,18);LAM=F(43,1000);NU=F(61,10000);C=LAM/NU;RHO=C*C;OMEGA=32;TAILV=1600
CASES=[{'name':'strong','D':F(1693,200),'floor':302,'scale':F(1),'target':F(17,50)}, {'name':'relaxed','D':F(8),'floor':283,'scale':F(3,2),'target':F(7,20)}]
def aa(q):
 q=F(q);return arb(q.numerator)/q.denominator

def pp(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def normalpeak(v):
 a=aa(ETA)*arb(v).sqrt();lo=F(-10);hi=F(10)
 def sg(z):
  z=aa(z);return (-z*z/2).exp()/(2*arb.pi()).sqrt()-a*(-z/arb(2).sqrt()).erfc()/2
 assert sg(lo)>0 and sg(hi)<0
 for _ in range(64):
  m=(lo+hi)/2;y=sg(m)
  if y>0:lo=m
  else:assert y<0;hi=m
 z=arb(aa((lo+hi)/2),aa((hi-lo)/2));h=(-(a+z)**2/2).exp()/((2*arb.pi()).sqrt()*a)
 return h.upper(),[str(lo),str(hi)]

def edgepeak(v,T):
 a=aa(ETA)*arb(v).sqrt();k=T/(6*arb(v)*arb(v).sqrt())
 assert k>0 and k<aa(F(1,10))
 def vals(z):
  z=aa(z) if isinstance(z,F) else z
  ph=(-z*z/2).exp()/(2*arb.pi()).sqrt()
  fg=(a*a/2-a*z).exp()*((a-z)/arb(2).sqrt()).erfc()/2
  P=z*z+a*z+a*a-1
  val=(1+k*a**3)*fg-k*P*ph
  der=(1+k*a**3)*(ph-a*fg)-k*(2*z+a-z*P)*ph
  return val,der
 assert vals(F(0))[0]>0 and vals(F(0))[1]>0 and vals(F(2))[1]<0
 lo=F(0);hi=F(2)
 for _ in range(64):
  mid=(lo+hi)/2;sg=vals(mid)[1]
  if sg>0:lo=mid
  else:assert sg<0;hi=mid
 z=arb(aa((lo+hi)/2),aa((hi-lo)/2));val,_=vals(z)
 return val.upper(),[str(lo),str(hi)],pp(k)

def setup():
 assert sha(TABLE)==TSHA
 table=json.loads(TABLE.read_text());ker=json.loads((HERE/'bandlimited_kernel_cells.json').read_text())
 cert=json.loads((HERE/'BANDWIDTH_KERNEL_INPUTS.json').read_text())
 assert sha(HERE/'bandlimited_kernel_cells.json')==cert['kernel_cells_sha256']
 assert ker['table_sha256']==TSHA and ker['omega']==OMEGA
 lc=json.loads((HERE/'LOCAL_FOURTH_CERTIFICATE.json').read_text())
 assert lc['status'].startswith('PASS_') and sha(HERE/'certify_local_fourth.py')==lc['source_sha256']
 cells=[]
 for c in ker['cells']:
  l=F(c['l']);u=F(c['u']);A=F(table['frequency_upper_values'][c['index']]);mass=aa(F(c['fourier_weight_upper']))
  cells.append((l,u,A,mass))
 assert len(cells)==15104 and cells[-1][1]==OMEGA
 delta=arb(OMEGA)/(2*arb.pi());excess=1/(delta*(1-(-aa(ETA)/delta).exp()))-1/aa(ETA)
 return cells,excess.upper()

def run(case,cells,excess):
 st=time.monotonic();D=case['D'];target=case['target'];low=[];near=[];highA=arb(0);farM=arb(0)
 for l,u,A,mass in cells:
  cap=aa(min(F(1),case['scale']*A))
  if u<=F(1,4):
   small=(RHO*u*u/2<1)
   Rfac=aa(u**4*(F(1,3)+C/(8*(1-RHO*u*u/2)))) if small else None
   low.append((aa(l*l/2),Rfac,aa(u**6/72),aa(u**3/6),cap,mass))
  else:
   highA+=mass*cap
   if l<1:near.append((aa(l*l/2),aa(u**3/6),mass))
   else:farM+=mass
 highA=highA.upper();farM=farM.upper();rows=[]
 for vl in range(case['floor'],TAILV):
  vu=vl+1;T=aa((LAM*vu-D)/NU);err=arb(0)
  for l2,Rfac,u6,u3,cap,mass in low:
   gauss=(-vl*l2).exp();refamp=gauss*(1+(T*u3)**2).sqrt();bound=cap+refamp
   if Rfac is not None:
    lowbound=gauss*((T*Rfac).expm1()+T*T*u6)
    bound=bound.upper().min(lowbound.upper())
   err+=mass*bound
  nearG=sum((mass*(-vl*l2).exp()*(1+(T*u3)**2).sqrt() for l2,u3,mass in near),arb(0))
  farG=farM*(-arb(vl)/2).exp()*(1+(T*OMEGA**3/6)**2).sqrt()
  h,hbr=normalpeak(vl);e,ebr,k= edgepeak(vl,T)
  refpeak=h.max(e)
  density=1/(2*arb.pi()*vl).sqrt()+T/(3*arb.pi()*vl*vl)
  baseline=refpeak+excess*density
  total=(baseline+err+highA+nearG+farG).upper()
  row={'V_interval':[vl,vu],'T_upper':pp(T),'gaussian_peak_upper':pp(h),'gaussian_root':hbr,'edgeworth_peak_upper':pp(e),'edgeworth_root':ebr,'dimensionless_cubic_coefficient':k,'reference_peak_upper':pp(refpeak),'area_error_upper':pp(excess*density),'low_error_upper':pp(err),'high_A_upper':pp(highA),'near_G_upper':pp(nearG),'far_G_upper':pp(farG),'total_upper':pp(total)}
  rows.append(row)
  if not total<aa(target):
   (HERE/('EDGEWORTH_FAILURE_'+case['name']+'.json')).write_text(json.dumps(row,indent=2)+'\n');raise AssertionError((case['name'],vl,str(total)))
  if len(rows)%200==0:print(json.dumps({'case':case['name'],'boxes':len(rows),'last_V':vl,'seconds':time.monotonic()-st}),flush=True)
 h,_=normalpeak(TAILV);tail=h+aa(F(28,25)/NU*(LAM*TAILV-D))/(arb(TAILV)**3).sqrt()
 assert TAILV>3*D/LAM and tail<aa(target)
 worst=max(rows,key=lambda r:F(r['total_upper']['upper']))
 return {'case':case['name'],'status':'PASS_CONTINUOUS_V_DOMAIN_WITH_SIGNED_REFERENCE','parameters':{k:(str(v) if isinstance(v,F) else v) for k,v in case.items()},'variance_boxes':len(rows),'finite_maximum':worst,'infinite_variance_tail':pp(tail),'rows':rows,'seconds':time.monotonic()-st}

def main():
 st=time.monotonic();cells,excess=setup();results=[run(case,cells,excess) for case in CASES]
 out={'status':'PASS_EDGEWORTH_CL_TRACE_CONDITIONAL_ON_EXPLICIT_LISTED_INPUTS','scope':'actual divisor probability measure; full threshold atom; no new p endpoint','parameters':{'eta':str(ETA),'lambda':str(LAM),'nu':str(NU),'c':str(C),'rho':str(RHO),'omega':OMEGA,'fourth_moment_per_abs_third':'8','V_tail':TAILV},'frequency_cells':len(cells),'cases':results,'bindings':{'prior_A_table':str(TABLE),'prior_A_sha256':TSHA,'local_fourth_sha256':sha(HERE/'LOCAL_FOURTH_CERTIFICATE.json'),'kernel_weights_sha256':sha(HERE/'bandlimited_kernel_cells.json'),'kernel_producer_sha256':sha(HERE/'kernel_cells.py'),'source_sha256':sha(__file__)},'environment':{'python':sys.version,'platform':platform.platform(),'arb_precision_bits':ctx.prec},'seconds':time.monotonic()-st}
 (HERE/'EDGEWORTH_TRACE_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n')
 for case in results:
  print(json.dumps({k:v for k,v in case.items() if k!='rows'},indent=2),flush=True)
if __name__=='__main__':main()
