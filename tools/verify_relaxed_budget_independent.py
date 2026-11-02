#!/usr/bin/env python3
"""Independent relaxed-budget replay: certified frequency rescaling, Mills brackets.
No agent producer imports. Four read-only CPU workers; shared Arb disclosed.
"""
from pathlib import Path
import os
from fractions import Fraction as F
import sys,math,json,hashlib,time,argparse,multiprocessing as mp
R=Path(os.environ['MARKOFF_WORKDIR']);S=R/'agents/gaussian_fixed_tail_certificate_round34';O=R/'root_logs'
sys.path.insert(0,str(R/'agents/attachment_audit/d58/deps'))
from flint import arb,acb,ctx
ctx.prec=256
ETA=F(1,18);BETA=F(40,9);EPS=F(63,1000);DELTA=F(829,9180)
LAM=F(43,1000);NU=F(61,10000);OLD_D=F(1693,200);D=F(8);P=10**153;XX=F(P,2);RELAX=F(3,2)
CC=LAM/NU;RHO=CC*CC;FLOOR=F(283);GAMMA=F(2,5);OMEGA=512;A=14000
G=R/'agents/selected_resource_moment_round38'
TABLE=DATA=GRIDS=BUDGET=None;LOW=NEAR=HIGHA=FARG=TAIL=None
def aa(x):
 x=F(x);return arb(x.numerator)/x.denominator
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def pack(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def floorroot(n,k):
 lo=0;hi=1<<((n.bit_length()+k-1)//k)
 while hi-lo>1:
  mid=(lo+hi)//2
  if mid**k<=n:lo=mid
  else:hi=mid
 assert lo**k<=n<(lo+1)**k
 return lo
def local_root(q):
 sc=1<<230;n=floorroot(sc**18//q,18)
 assert q*n**18<=sc**18<q*(n+1)**18
 return arb(aa(F(2*n+1,2*sc)),aa(F(1,2*sc)))
def frequency_grid():
 return ([(F(k,16384),F(k+1,16384)) for k in range(4096)]
  +[(F(k,1024),F(k+1,1024)) for k in range(256,8192)]
  +[(F(k,128),F(k+1,128)) for k in range(1024,65536)])
def variance_grid():
 return [(F(k),F(k+1)) for k in range(283,1600)]
def setup():
 global TABLE,DATA,GRIDS,BUDGET
 TABLE=json.loads((S/'EULER_FREQUENCY_CERTIFICATE.json').read_text());params=TABLE['parameters']
 for key,want in {'eta':ETA,'beta':BETA,'epsilon':EPS,'delta':DELTA,'lambda':LAM,'nu':NU,'D':OLD_D,'A':A,'X':XX,'Omega':OMEGA}.items():
  assert F(params[key])==want,(key,params[key],want)
 GRIDS=frequency_grid();assert len(GRIDS)==len(TABLE['frequency_upper_values'])==76544
 BUDGET=arb(A).log()+aa(DELTA-EPS)*aa(XX).log()
 assert BUDGET-arb(43540).log()>aa(OLD_D)
 assert BUDGET-arb(43540).log()-aa(RELAX).log()>aa(D)
 assert 1+F(41,100)+F(41,100)**2/2+F(41,100)**3/6>RELAX
 assert OLD_D-F(41,100)>D and RELAX>=1
 ps=[q for q in range(2,4099) if all(q%d for d in range(2,math.isqrt(q)+1))]
 assert len(ps)==564 and ps==[x['q'] for x in TABLE['local_tails']]
 DATA=[];count=0
 for rec in TABLE['local_tails']:
  q=rec['q'];a0=rec['first_omitted_exponent'];h=arb(q).log();r=local_root(q);c=arb(q-1)/q
  assert (-aa(EPS)*a0*h).exp()*(1+c*r/(1-r))<1
  vals=[]
  for a in range(1,a0):
   weight=(-aa(EPS)*a*h).exp()
   sr=r*(1-(a+1)*r**a+a*r**(a+1))/(1-r)**2
   assert sr>0
   vals.append((a,weight,weight*c*h*sr))
  DATA.append((h,r,c,vals));count+=len(vals)
 assert count==1355 and (-aa(EPS)*arb(4099).log()).exp()*(1+local_root(4099))<1
 resource=json.loads((R/'control/SELECTED_RESOURCE_MOMENT_ROOT_ACCEPTANCE.json').read_text())
 assert resource['status']=='ACCEPTED_GENERAL_MOMENT_RESOURCE_INEQUALITY_NO_NEW_P_ENDPOINT'
 assert F(resource['C_upper_independent'])<164 and (D-NU*164)/(LAM-3*NU)>FLOOR
def frequency_worker(chunk):
 begin,end=chunk;fresh=trivial=0;ratio=arb(0);start=time.monotonic()
 invbudget=(-BUDGET).exp()
 for i in range(begin,end):
  saved=F(TABLE['frequency_upper_values'][i]);assert 0<saved<=1
  if saved==1:trivial+=1;continue
  lo,hi=GRIDS[i];w=aa((lo+hi)/2);rad=aa((hi-lo)/2);prod=arb(1)
  for h,r,c,vals in DATA:
   z=acb(0,h*w).exp()*r;den=1-z;assert den.abs_lower()>0
   best=arb(1)
   for a,wt,der in vals:
    v=1+c*z*(1-z**a)/den
    upper=(wt*v.abs_upper()+der*rad).upper()
    best=best.max(upper).upper()
   prod*=best
  bound=(prod*invbudget).upper()
  assert bound<=aa(saved),(i,str(bound),str(saved))
  ratio=ratio.max((bound/aa(saved)).upper()).upper();fresh+=1
 return {'begin':begin,'end':end,'fresh_nontrivial_cells':fresh,'trivial_unit_cells':trivial,'largest_ratio_upper':str(ratio.fmpq()),'seconds':time.monotonic()-start}
def mills_ratio(u):
 # Taylor bounds on e^(-t^2/2), integrated against e^(-u*t).
 assert u>50
 inv=1/u;inv2=inv*inv;term=inv;upper=term
 for k in range(1,21):
  term*=-(2*k-1)*inv2;upper+=term
 next_term=term*(-41)*inv2
 lower=upper+next_term
 assert lower.lower()>0 and lower.lower()<=upper.upper()
 lo=lower.lower();hi=upper.upper()
 return arb((lo+hi)/2,((hi-lo)/2).upper())
def normal_peak(v):
 sig=aa(v).sqrt();rt2=arb(2).sqrt();rtpi=(2*arb.pi()).sqrt()
 def parts(t):
  t=t if isinstance(t,arb) else aa(t)
  w=(t-aa(ETA)*aa(v))/sig
  first=(aa(ETA)**2*aa(v)/2-aa(ETA)*t).exp()*(1+(w/rt2).erf())/2
  u=(t+aa(BETA)*aa(v))/sig
  second=(-t*t/(2*aa(v))).exp()/rtpi*mills_ratio(u)
  return first,second
 def der(t):
  x,y=parts(t);return aa(BETA)*y-aa(ETA)*x
 lo=F(-64);hi=F(64);assert der(lo)>0 and der(hi)<0
 for _ in range(60):
  mid=(lo+hi)/2;d=der(mid)
  if d>0:lo=mid
  else:assert d<0;hi=mid
 x,y=parts(arb(aa((lo+hi)/2),aa((hi-lo)/2)))
 return (x+y).upper()
def half_peak(v):
 a=aa(v).sqrt()/18;rt2=arb(2).sqrt();rtpi=(2*arb.pi()).sqrt()
 def der(w):
  w=aa(w);return (-w*w/2).exp()/rtpi-a*(1+(w/rt2).erf())/2
 lo=F(-10);hi=F(10);assert der(lo)>0 and der(hi)<0
 for _ in range(60):
  mid=(lo+hi)/2
  if der(mid)>0:lo=mid
  else:assert der(mid)<0;hi=mid
 w=arb(aa((lo+hi)/2),aa((hi-lo)/2))
 return ((-a*a/2-a*w).exp()*(1+(w/rt2).erf())/2).upper()
def variance_setup():
 global LOW,NEAR,HIGHA,FARG,TAIL
 LOW=[];NEAR=[];high=arb(0);pi=arb.pi()
 for (lo,hi),sv in zip(GRIDS,TABLE['frequency_upper_values']):
  a=aa(min(F(1),RELAX*F(sv)));l=aa(lo);u=aa(hi)
  mass=(aa(hi-lo)*aa(ETA+BETA)/(pi*(aa(ETA*ETA+lo*lo)*aa(BETA*BETA+lo*lo)).sqrt())).upper()
  if hi<=F(1,4):LOW.append((l*l/2,u**3/3,aa(RHO)*u*u/2,u**3/6+aa(CC)*u**4/8,a,mass))
  else:
   high+=mass*a
   if lo<1:NEAR.append((l*l/2,mass))
 HIGHA=high.upper()
 # Independent analytic tail for Gaussian beyond1, rather than producer finite mass sum.
 FARG=(aa(ETA+BETA)/pi*(-aa(FLOOR)/2).exp()).upper()
 TAIL=(aa(ETA+BETA)/(OMEGA*pi)*(1+(-aa(FLOOR)*OMEGA**2/2).exp())).upper()
 assert len(LOW)==4096 and len(NEAR)==768
def variance_worker(chunk):
 begin,end=chunk;grid=variance_grid();rows=[];maximum=arb(0);worst=None
 critical=3*D/LAM
 for idx in range(begin,end):
  vl,vu=grid[idx];v=aa(vl);t=aa((LAM*vu-D)/NU);low=arb(0)
  for l2,u3,rhou2,pref,a,mass in LOW:
   ev=-v*l2
   err=t*pref*(ev+t*u3+rhou2).exp()
   cap=a+ev.exp()
   low+=mass*err.upper().min(cap.upper()).upper()
  near=sum((mass*(-v*l2).exp() for l2,mass in NEAR),arb(0))
  rep=normal_peak(vl)+low+HIGHA+near+FARG+TAIL
  arg=min(vu,max(vl,critical))
  be=half_peak(vl)+aa(F(28,25)/NU)*aa(LAM*arg-D)/(aa(arg)**3).sqrt()
  total=rep.upper().min(be.upper()).upper()
  assert total<aa(GAMMA),(vl,vu,str(total))
  if total>maximum:maximum=total;worst=[str(vl),str(vu)]
  rows.append({'V_interval':[str(vl),str(vu)],'replacement_upper':str(rep.upper().fmpq()),'BE_upper':str(be.upper().fmpq()),'total_upper':str(total.fmpq())})
 return {'begin':begin,'end':end,'maximum_upper':str(maximum.fmpq()),'maximum_box':worst,'rows':rows}
def chunks(n,workers):return [(i*n//workers,(i+1)*n//workers) for i in range(workers)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
 assert 1<=args.workers<=8;start=time.monotonic();setup()
 h=sha(S/'EULER_FREQUENCY_CERTIFICATE.json')
 assert h=='0bfce18d9fbe285fd4bd4334c5c0d477b954d5206497a2ac392119f6f8a4c656'
 freq_receipt=O/'RESOURCE_GAUSSIAN_FREQUENCY_ROOT_REVIEW.json'
 prior=json.loads(freq_receipt.read_text())
 assert prior['table_sha256']==h and prior['status']=='PASS_INDEPENDENT_FIXED_FREQUENCY_ALL76544_CELLS'
 accept=R/'control/RESOURCE_GAUSSIAN_SIMPLIFICATION_ROOT_ACCEPTANCE.json'
 accepted=json.loads(accept.read_text())
 assert accepted['status']=='ACCEPTED_SIMPLIFIED_INTERNAL_PROOF_CHAIN_SAME_P_ENDPOINT'
 assert accepted['gamma_proved']=='177/500'
 assert (D-164*NU)/(LAM-3*NU)==F(69996,247)>FLOOR
 assert F(0)<DELTA-EPS<F(1,4) and F(4,5)>F(14,15)**4
 variance_setup();grid=variance_grid();n=len(grid);assert n==1317
 with mp.get_context('fork').Pool(args.workers) as pool:res=pool.map(variance_worker,chunks(n,args.workers))
 rows=[row for x in res for row in x['rows']]
 assert len(rows)==n and [x['V_interval'] for x in rows]==[[str(a),str(b)] for a,b in grid]
 maximum=max(F(x['total_upper']) for x in rows)
 repmax=max(F(x['replacement_upper']) for x in rows)
 which=next(x['V_interval'] for x in rows if F(x['total_upper'])==maximum)
 assert maximum<GAMMA and 1600>3*D/LAM
 infinity=half_peak(1600)+aa(F(28,25)/NU)*aa(LAM*1600-D)/(arb(1600)**3).sqrt()
 assert infinity<aa(GAMMA)
 precheck=O/'BUDGET_RELAXATION_PRECHECK_AND_GENERAL_LEMMA.md'
 out={'status':'PASS_INDEPENDENT_RELAXED_BUDGET_ALL1317_VARIANCE_BOXES_AND_TAILS',
      'no_new_p_acceptance':True,'gamma':str(GAMMA),'D':str(D),'variance_floor':283,'variance_boxes':n,
      'frequency_relaxation_factor':str(RELAX),'frequency_rule':'min(1,(3/2)*old_certified_upper)',
      'no_new_frequency_replay_claimed':True,'rho':str(RHO),'fourth_order_c':str(CC),
      'maximum_finite_upper':str(maximum),'maximum_finite_display':str(aa(maximum)),'maximum_box':which,
      'maximum_replacement_upper':str(repmax),'all_replacement_only_bounds_below_gamma':repmax<GAMMA,
      'infinite_variance_tail_upper':pack(infinity),'frequency_tail_upper':pack(TAIL),
      'high_A_upper':pack(HIGHA),'Gaussian_far_upper':pack(FARG),'rows':rows,
      'frequency_receipt':str(freq_receipt.relative_to(R)),'frequency_receipt_sha256':sha(freq_receipt),
      'table_sha256':h,'parent_acceptance_sha256':sha(accept),'precheck_sha256':sha(precheck),
      'verifier_sha256':sha(__file__),'seconds':time.monotonic()-start,
      'method':'Root prior independent original-CDF/Mills-Taylor implementation; fresh60-step roots and full1317-box rectangle integration after rigorous rational cap rescaling. Independent of new agent producer, shared Arb256.'}
 (O/'RELAXED_BUDGET_VARIANCE_ROOT_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2),flush=True)
 assert sha(S/'EULER_FREQUENCY_CERTIFICATE.json')==h
if __name__=='__main__':main()
