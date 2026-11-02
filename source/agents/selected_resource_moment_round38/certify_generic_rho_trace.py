#!/usr/bin/env python3
"""Fixed simplification: generic rho, new resource floor, old fixed frequency table.
No producer/helper imports. This does not independently recalculate the frequency table.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,math,json,hashlib,time,argparse,multiprocessing as mp
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb,ctx
ctx.prec=256
ETA=F(1,18);BETA=F(40,9);EPS=F(63,1000);DELTA=F(829,9180)
LAM=F(43,1000);NU=F(61,10000);D=F(1693,200);C=LAM/NU;RHO=C*C
FLOOR=302;END=1600;GAMMA=F(189,500);OMEGA=512
TABLE_PATH=ROOT/'agents/gaussian_fixed_tail_certificate_round34/EULER_FREQUENCY_CERTIFICATE.json'
TABLE_SHA='0bfce18d9fbe285fd4bd4334c5c0d477b954d5206497a2ac392119f6f8a4c656'
LOW=NEAR=HIGH_A=FAR_G=TAIL=None
def aa(q):q=F(q);return arb(q.numerator)/q.denominator
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def pp(x):return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def grid():
 return ([(F(k,16384),F(k+1,16384)) for k in range(4096)]
 +[(F(k,1024),F(k+1,1024)) for k in range(256,8192)]
 +[(F(k,128),F(k+1,128)) for k in range(1024,65536)])
def setup():
 global LOW,NEAR,HIGH_A,FAR_G,TAIL
 assert sha(TABLE_PATH)==TABLE_SHA
 data=json.loads(TABLE_PATH.read_text());pa=data['parameters']
 for key,want in {'eta':ETA,'beta':BETA,'epsilon':EPS,'delta':DELTA,'A':14000,'X':F(10**153,2),'Omega':512}.items():assert F(pa[key])==want
 gs=grid();assert len(gs)==len(data['frequency_upper_values'])==76544
 rc=json.loads((HERE/'RESOURCE_CERTIFICATE.json').read_text());assert F(rc['C_upper'])<164 and rc['safe_V_floor']==302
 assert (D-NU*164)/(LAM-3*NU)>FLOOR
 assert arb(14000).log()+aa(DELTA-EPS)*aa(F(10**153,2)).log()-arb(43540).log()>aa(D)
 LOW=[];NEAR=[];hs=arb(0);pi=arb.pi()
 for (l,u),sv in zip(gs,data['frequency_upper_values']):
  av=F(sv);assert 0<av<=1;a=aa(av)
  mass=(aa(u-l)*aa(ETA+BETA)/(pi*(aa(ETA*ETA+l*l)*aa(BETA*BETA+l*l)).sqrt())).upper()
  if u<=F(1,4):LOW.append((aa(l*l/2),aa(u**3/3),aa(RHO*u*u/2),aa(u**3/6+C*u**4/8),a,mass))
  else:
   hs+=mass*a
   if l<1:NEAR.append((aa(l*l/2),mass))
 HIGH_A=hs.upper()
 FAR_G=(aa(ETA+BETA)/pi*(-arb(FLOOR)/2).exp()).upper()
 TAIL=(aa(ETA+BETA)/(OMEGA*pi)*(1+(-arb(FLOOR)*OMEGA**2/2).exp())).upper()
 assert len(LOW)==4096 and len(NEAR)==768
def peak(v):
 sig=(2*arb(v)).sqrt()
 def parts(L):
  x=L if isinstance(L,arb) else aa(L)
  p=(aa(ETA*ETA*v/2)-aa(ETA)*x).exp()*((aa(ETA*v)-x)/sig).erfc()/2
  q=(aa(BETA*BETA*v/2)+aa(BETA)*x).exp()*((aa(BETA*v)+x)/sig).erfc()/2
  return p,q
 def der(L):p,q=parts(L);return aa(BETA)*q-aa(ETA)*p
 lo=F(-100);hi=F(100);assert der(lo)>0 and der(hi)<0
 for _ in range(48):
  mid=(lo+hi)/2;value=der(mid)
  if value>0:lo=mid
  else:assert value<0;hi=mid
 p,q=parts(arb(aa((lo+hi)/2),aa((hi-lo)/2)))
 return (p+q).upper(),[str(lo),str(hi)]
def worker(task):
 begin,end=task;rows=[];st=time.monotonic()
 for vl in range(begin,end):
  vu=vl+1;t=aa((LAM*vu-D)/NU);low=arb(0)
  for l2,u3,rhou2,pref,av,mass in LOW:
   ev=-vl*l2
   err=t*pref*(ev+t*u3+rhou2).exp()
   cap=av+ev.exp()
   low+=mass*err.upper().min(cap.upper()).upper()
  near=sum((mass*(-vl*l2).exp() for l2,mass in NEAR),arb(0))
  peak_bound,bracket=peak(vl)
  full=(peak_bound+low+HIGH_A+near+FAR_G+TAIL).upper()
  row={'V_interval':[str(vl),str(vu)],'T_upper':pp(t),'Gaussian_root_interval':bracket,'Gaussian_peak_upper':pp(peak_bound),'low_error_upper':pp(low.upper()),'near_G_upper':pp(near.upper()),'total_upper':pp(full)}
  if not full<aa(GAMMA):
   (HERE/f'FAILURE_{vl}.json').write_text(json.dumps(row,indent=2)+'\n');raise AssertionError((vl,str(full)))
  rows.append(row)
  if (vl-begin+1)%100==0:print(json.dumps({'worker_interval':[begin,end],'completed':vl-begin+1,'latest_V':vl,'seconds':time.monotonic()-st}),flush=True)
 return rows
def main():
 st=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);args=ap.parse_args();assert 1<=args.workers<=4
 setup();tasks=[(FLOOR+i*(END-FLOOR)//args.workers,FLOOR+(i+1)*(END-FLOOR)//args.workers) for i in range(args.workers)]
 with mp.get_context('fork').Pool(args.workers) as pool:rows=[r for part in pool.map(worker,tasks) for r in part]
 assert [r['V_interval'] for r in rows]==[[str(v),str(v+1)] for v in range(FLOOR,END)]
 maxrow=max(rows,key=lambda x:F(x['total_upper']['upper']));finite=F(maxrow['total_upper']['upper'])
 tail=(18/(2*arb.pi()).sqrt()+aa(F(28,25)*LAM/NU))/arb(END).sqrt();assert tail<aa(GAMMA)
 global_max=max(finite,F(str(tail.upper().fmpq())))
 out={'status':'GENERIC_RHO_RESOURCE_TRACE_CERTIFIED_PENDING_ROOT_ACCEPTANCE','scope':'Same fixed frequency table; generic rho=(lambda/nu)^2; new independent resource floor302. No giant selected-variance table dependency. Frequency independent replay remains root responsibility.','parameters':{'eta':str(ETA),'beta':str(BETA),'epsilon':str(EPS),'delta':str(DELTA),'lambda':str(LAM),'nu':str(NU),'D':str(D),'rho':str(RHO),'fourth_order_coefficient_c':str(C),'A':14000,'X':str(F(10**153,2)),'Omega':OMEGA,'gamma':str(GAMMA),'V_floor':FLOOR,'V_tail':END},'variance_boxes':len(rows),'frequency_cells':76544,'finite_maximum':maxrow,'global_upper':str(global_max),'global_display_only':str(aa(global_max)),'high_A_upper':pp(HIGH_A),'far_G_upper':pp(FAR_G),'frequency_tail_upper':pp(TAIL),'infinite_variance_tail_upper':pp(tail),'rows':rows,'bindings':{str(TABLE_PATH.relative_to(ROOT)):TABLE_SHA,str((HERE/'RESOURCE_CERTIFICATE.json').relative_to(ROOT)):sha(HERE/'RESOURCE_CERTIFICATE.json')},'producer_sha256':sha(__file__),'seconds':time.monotonic()-st,'external_theorem':'Inherited independent nonidentical Berry-Esseen C=.56 only for V>=1600 tail.','arithmetic_backend':'Arb256; no predecessor code/helper imports.'}
 (HERE/'GENERIC_RHO_TRACE_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['rows','finite_maximum']},indent=2),flush=True);print('finite maximum',maxrow['V_interval'],maxrow['total_upper']['display_only'],flush=True)
 assert sha(TABLE_PATH)==TABLE_SHA
if __name__=='__main__':main()
