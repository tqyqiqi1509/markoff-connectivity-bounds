#!/usr/bin/env python3
"""Independent resource-floor replay: closed sums, dyadic roots, Mills brackets.
No agent producer imports. Four read-only CPU workers; shared Arb disclosed.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,math,json,hashlib,time,argparse,multiprocessing as mp
R=Path(__file__).resolve().parents[1];S=R/'agents/gaussian_fixed_tail_certificate_round34';O=R/'root_logs'
sys.path.insert(0,str(R/'agents/attachment_audit/d58/deps'))
from flint import arb,acb,ctx
ctx.prec=256
ETA=F(1,18);BETA=F(40,9);EPS=F(63,1000);DELTA=F(829,9180)
LAM=F(43,1000);NU=F(61,10000);D=F(1693,200);P=10**153;XX=F(P,2)
CC=LAM/NU;RHO=CC*CC;FLOOR=F(302);GAMMA=F(177,500);OMEGA=512;A=14000
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
 return [(F(k),F(k+1)) for k in range(302,1600)]
def setup():
 global TABLE,DATA,GRIDS,BUDGET
 TABLE=json.loads((S/'EULER_FREQUENCY_CERTIFICATE.json').read_text());params=TABLE['parameters']
 for key,want in {'eta':ETA,'beta':BETA,'epsilon':EPS,'delta':DELTA,'lambda':LAM,'nu':NU,'D':D,'A':A,'X':XX,'Omega':OMEGA}.items():
  assert F(params[key])==want,(key,params[key],want)
 GRIDS=frequency_grid();assert len(GRIDS)==len(TABLE['frequency_upper_values'])==76544
 BUDGET=arb(A).log()+aa(DELTA-EPS)*aa(XX).log()
 assert BUDGET-arb(43540).log()>aa(D)
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
  a=aa(F(sv));l=aa(lo);u=aa(hi)
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
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['frequency','variance','all'],default='all');ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
 assert 1<=args.workers<=8;start=time.monotonic();setup();h=sha(S/'EULER_FREQUENCY_CERTIFICATE.json')
 if args.stage in ['frequency','all']:
  with mp.get_context('fork').Pool(args.workers) as pool:res=pool.map(frequency_worker,chunks(76544,args.workers))
  assert sum(x['fresh_nontrivial_cells']+x['trivial_unit_cells'] for x in res)==76544
  out={'status':'PASS_INDEPENDENT_FIXED_FREQUENCY_ALL76544_CELLS','table_sha256':h,'local_prime_count':564,'finite_exponents':1355,'workers':res,'verifier_sha256':sha(__file__),'seconds':time.monotonic()-start,'method':'Direct complex exponential at every midpoint, closed finite geometric sum and differentiated closed sum, exact230-bit dyadic18th roots; full prime/exponent tails. Shared Arb256.'}
  (O/'RESOURCE_GAUSSIAN_FREQUENCY_ROOT_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
 if args.stage in ['variance','all']:
  freq_receipt=O/'RESOURCE_GAUSSIAN_FREQUENCY_ROOT_REVIEW.json'
  if not freq_receipt.exists():freq_receipt=O/'FIXED_GAUSSIAN_FREQUENCY_ROOT_REVIEW.json'
  prior=json.loads(freq_receipt.read_text());assert prior['table_sha256']==h
  assert prior['status']=='PASS_INDEPENDENT_FIXED_FREQUENCY_ALL76544_CELLS'
  cert=json.loads((G/'GENERIC_RHO_TRACE_CERTIFICATE.json').read_text());assert cert['bindings']['agents/gaussian_fixed_tail_certificate_round34/EULER_FREQUENCY_CERTIFICATE.json']==h
  assert [x['V_interval'] for x in cert['rows']]==[[str(a),str(b)] for a,b in variance_grid()]
  variance_setup()
  with mp.get_context('fork').Pool(args.workers) as pool:res=pool.map(variance_worker,chunks(1298,args.workers))
  rows=[row for x in res for row in x['rows']];assert len(rows)==1298
  maximum=max(F(x['total_upper']) for x in rows);which=next(x['V_interval'] for x in rows if F(x['total_upper'])==maximum)
  assert 1600>3*D/LAM
  infinity=half_peak(1600)+aa(F(28,25)/NU)*aa(LAM*1600-D)/(arb(1600)**3).sqrt();assert infinity<aa(GAMMA)
  out={'status':'PASS_INDEPENDENT_RESOURCE_GAUSSIAN_ALL1298_VARIANCE_BOXES_AND_SHARP_TAIL','gamma':'177/500','variance_floor':302,'variance_boxes':1298,'rho':str(RHO),'fourth_order_c':str(CC),'maximum_finite_upper':str(maximum),'maximum_finite_display':str(aa(maximum)),'maximum_box':which,'infinite_variance_tail_upper':pack(infinity),'frequency_tail_upper':pack(TAIL),'high_A_upper':pack(HIGHA),'Gaussian_far_upper':pack(FARG),'rows':rows,'frequency_receipt':str(freq_receipt.relative_to(R)),'frequency_receipt_sha256':sha(freq_receipt),'table_sha256':h,'certificate_sha256':sha(G/'GENERIC_RHO_TRACE_CERTIFICATE.json'),'verifier_sha256':sha(__file__),'seconds':time.monotonic()-start,'method':'Original CDF first term plus rigorous alternating integrated-Taylor Mills ratio for second term; fresh60-step roots, full rectangle integration, independent analytic high-G tail. Shared Arb256.'}
  (O/'RESOURCE_GAUSSIAN_VARIANCE_ROOT_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2),flush=True)
 assert sha(S/'EULER_FREQUENCY_CERTIFICATE.json')==h
if __name__=='__main__':main()
