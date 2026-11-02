#!/usr/bin/env python3
"""Arb/exact producer for selected real moments and a uniform trace endpoint.

Reuses existing elementary integer/Fraction helpers only for old gates and
endpoint roots. It is not an independent replay of those helpers.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,importlib.util,argparse
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb,ctx
ctx.prec=192
PRIOR=ROOT/'agents/analytic_window_round4'
spec=importlib.util.spec_from_file_location('existing_elementary_trace',PRIOR/'certify_trace_concentration.py')
OLD=importlib.util.module_from_spec(spec);spec.loader.exec_module(OLD)
HLP=OLD.HLP;pack=OLD.pack
P0=3*10**153;ETA=F(1,18);DELTA=F(829,9180);DP=F(63,1000)
TT=F(43,1000);UU=F(61,10000);DD=F(8747,1000);BE=F(14,25)

def aa(q):
 q=F(q);return arb(q.numerator)/q.denominator
def rat(q):return F(int(q['numerator']),int(q['denominator']))
def arpack(x):
 assert x.is_finite()
 return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def euler():
 rows=[];product=arb(1);t=aa(TT);u=aa(UU);count=0
 for q in HLP.primes_before(4099):
  lq=arb(q).log();r=(-aa(ETA)*lq).exp();w=(-aa(DP)*lq).exp();c=arb(q-1)/q
  binf=1+c*r/(1-r);wp=arb(1);rp=arb(1);weights=[arb(1)];values=[];local=arb(1);amax=0
  a=0
  while True:
   a+=1;wp*=w;rp*=r
   if wp*binf<1:break
   weights.append(c*rp);B=sum(weights,arb(0));prob=[v/B for v in weights]
   xs=[j*lq for j in range(a+1)];mu=sum((pr*x for pr,x in zip(prob,xs)),arb(0))
   var=sum((pr*(x-mu)**2 for pr,x in zip(prob,xs)),arb(0))
   third=sum((pr*abs(x-mu)**3 for pr,x in zip(prob,xs)),arb(0))
   W=t*var-u*third;Wlo=W.lower() if W>0 else arb(0)
   value=wp*B*(-Wlo).exp();upper=value.upper()
   if upper>local:local=upper;amax=a
   values.append({'a':a,'B':arpack(B),'mean':arpack(mu),'variance':arpack(var),
                  'third_absolute_centered':arpack(third),'W_lower':arpack(Wlo),
                  'local_value_upper':arpack(upper)})
   count+=1;assert a<1000
  product*=local
  rows.append({'q':q,'maximizer_of_enclosures':amax,'local_upper':arpack(local),'values':values,
               'tail_first_exponent':a,'tail_guard':arpack(wp*binf)})
 lq=arb(4099).log();prime_tail=(-aa(DP)*lq).exp()*(1+(-aa(ETA)*lq).exp())
 assert prime_tail<1 and product<43540
 size_D=arb(18000).log()+aa(DELTA-DP)*(arb(P0)/2).log()-arb(43540).log()
 assert size_D>aa(DD)
 return {'status':'ARBITRARY_PRIME_EXPONENT_EULER_PASSED','precision_bits':ctx.prec,
         'eta':str(ETA),'delta_reference':str(DELTA),'delta_prime':str(DP),
         't':str(TT),'u':str(UU),'K_integer_cap':43540,'K_upper':arpack(product),
         'A':18000,'P0':str(P0),'D_lower':str(DD),'D_size_expression':arpack(size_D),
         'prime_tail_from':4099,'prime_tail_guard':arpack(prime_tail),
         'prime_count':len(rows),'finite_nonzero_exponent_count':count,'rows':rows}

def gaussian_max(v):
 a=aa(v).sqrt()/18;pi=arb.pi();sqrt2=arb(2).sqrt()
 def ratio(w):
  z=aa(w);return (-z*z/2).exp()/(2*pi).sqrt()/((1+(z/sqrt2).erf())/2)
 lo=F(-10);hi=F(10)
 assert ratio(lo)>a and ratio(hi)<a
 for _ in range(50):
  mid=(lo+hi)/2;rr=ratio(mid)
  if rr>a:lo=mid
  elif rr<a:hi=mid
  else:raise ArithmeticError('Increase precision: unseparated normal root sign')
 assert ratio(lo)>a and ratio(hi)<a
 w=arb(aa((lo+hi)/2),aa((hi-lo)/2).upper())
 H=(-(a+w)**2/2).exp()/((2*pi).sqrt()*a)
 return H,lo,hi

def variance_domain():
 assert DD/TT>203
 rows=[];maximum=arb(0);which=None;critical=3*DD/TT
 for vl in range(203,1200):
  vu=vl+1;H,lo,hi=gaussian_max(vl)
  vm=max(F(vl),min(F(vu),critical))
  er=2*aa(BE)/aa(UU)*(aa(TT)*aa(vm)-aa(DD))/(aa(vm)**3).sqrt()
  if er<0:er=arb(0)
  bound=H.upper()+er.upper()
  assert bound<aa(F(97,200))
  if bound>maximum:maximum=bound;which=vl
  rows.append({'variance_interval':[vl,vu],'normal_root_interval':[str(lo),str(hi)],
               'gaussian_H_upper':arpack(H.upper()),'error_argmax':str(vm),
               'BE_error_upper':arpack(er.upper()),'total_upper':arpack(bound)})
 tail=(18/(2*arb.pi()).sqrt()+2*aa(BE)*aa(TT)/aa(UU))/arb(1200).sqrt()
 assert tail<aa(F(97,200))
 gamma=F(97,200);sum_ratio=gamma+(1-gamma)*F(159842,18000**2+159842)
 assert sum_ratio<F(243,500)
 return {'status':'CONTINUOUS_ALL_VARIANCES_PASSED','variance_floor':str(DD/TT),
         'variance_boxes':len(rows),'maximum_box':which,'maximum_upper':arpack(maximum),
         'tail_from':1200,'tail_upper':arpack(tail),'gamma_cap':'97/200',
         'sum_ratio_upper':pack(sum_ratio),'sum_ratio_cap':'243/500','rows':rows}

def old_probability_gate():
 source=PRIOR/'SIZE_TRACE_CERTIFICATE.json';c=json.loads(source.read_text())
 for rec in c['dependencies']:
  p=Path(rec['path']);assert sha(p)==rec['sha256'] and p.stat().st_size==rec['bytes']
 e=c['euler'];assert e['integer_cap']==40300 and rat(e['K_upper'])<40300
 assert 40300**18360*2**281<183**18360*P0**281
 a,b=OLD.atan_bounds(F(1,5));c,d=OLD.atan_bounds(F(1,239));pi_l=16*a-4*d
 log_l=5*OLD.log_lower(F(2))+OLD.log_lower(F(10000,32*183))
 q=F(1411,10000);bound=F(1,4)/(pi_l*log_l);assert bound<q*q
 s=OLD.exp_minus_upper(F(5,24));k=7;r=1-k*q;assert k*q<=1<(k+1)*q
 gamma=q*sum((s**j for j in range(k)),F(0))+r*s**k
 ratio=gamma+(1-gamma)*F(159842,10000**2+159842);assert ratio<F(29,50)
 return {'status':'OLD_LP_REPROVED_AT_NEW_SCALE','prior_size_certificate_sha256':sha(source),
         'P0':str(P0),'effective_cap':183,'q_cap':str(q),'remainder':str(r),
         'q_squared_upper':pack(bound),'gamma_upper':pack(gamma),'sum_ratio_upper':pack(ratio),
         'sum_ratio_cap':'29/50'}

def endpoint():
 P=P0;H=F(8840);L=F(301);Cm=F(34195);CE=F(28924);K=F(159842);KL=F(1826);cc=F(11,10);C=30130
 ir=HLP.integer_root((P+1)**5,44);T=ir+int(ir**44!=(P+1)**5)
 assert (T-1)**44<(P+1)**5<=T**44
 safe=36864*2541933**4*2541935**2*T**6;cage=256*2541933**4*T**4
 assert safe<P and cage<P
 assert 21718**201960*2541935**660<int(H)**201960*2**119*P**541
 other=(729*(2*181265+3))**900*27**420*14556**6660*(P0+1)**673;other_rhs=(16*P0)**900
 assert other<other_rhs
 eL,eU=HLP.power_box(CE/H,15,2);mL,mU=HLP.power_box(Cm/H,15,2)
 xh=K/KL;yl=L/cc;xl=K/yl
 assert 0<xh<xl<eL<=eU<10000<18000<mL<=mU<K
 _,rp=HLP.power_box(K,-2,15)
 pts=[('left_limit',F(0),KL,H+L,H+L),('horizontal_to_hyperbola',xh,KL,H+L,H+L),
      ('low_count_crossover',xl,yl,H+L,H+L),
      ('high_mass_crossover',eU,K/eL,H+cc*K/eL,H+cc*K/eL),
      ('high_count_crossover',mU,K/mL,H+cc*K/mL,H*CE/Cm+cc*K/mL),
      ('right_edge',K,F(1),Cm*rp+cc,CE*rp+cc)]
 def split(A):
  A=F(A);_,ap=HLP.power_box(A,-2,15)
  return ('split_x_'+str(A),A,K/A,min(H,Cm*ap)+cc*K/A,min(H,CE*ap)+cc*K/A)
 s1=split(10000);s2=split(18000);rows=[]
 for name,items,factor in [('low',pts[:4]+[s1],F(1)),('middle',[s1,s2],F(29,50)),('high',[s2]+pts[4:],F(243,500))]:
  for label,x,y,U,Vmass in items:
   base=U**1510*Vmass**377*((x+y)*factor)**252;assert base<C**1887
   _,hi=HLP.power_box(base,1,1887)
   rows.append({'region':name,'point':label,'x_upper':pack(x),'y_upper':pack(y),
                'U_upper':pack(U),'V_mass_upper':pack(Vmass),'trace_ratio_used':str(factor),
                'F_upper':pack(hi),'F_display_only':float(hi)})
 left=(3*C)**142783*(P+1)**14616;right=2**81127*P**19068;assert left<right
 return {'status':'EXACT_MARKOFF_ENDPOINT_PASSED_PENDING_ROOT_REVIEW','P':str(P),'C':C,'ordinary_high_cap':int(H),
         'gates':{'safe8':True,'cage':True,'ordinary_transfer':True,'other63_at_3e153':True,'three_region_full_domain':True,'endpoint_integer_comparison':True},
         'other63_comparison_bits':[other.bit_length(),other_rhs.bit_length()],
         'endpoint_comparison_bits':[left.bit_length(),right.bit_length()],
         'point_bounds':rows,'prerequisites':{'T':str(T),'safe8_left':str(safe),'cage_left':str(cage)}}

def main():
 pa=argparse.ArgumentParser();pa.add_argument('--arithmetic-only',action='store_true');args=pa.parse_args()
 start=time.monotonic();out={'status':'PENDING'}
 out['euler']=euler();print('Euler PASS',out['euler']['prime_count'],out['euler']['finite_nonzero_exponent_count'],flush=True)
 out['variance_domain']=variance_domain();print('All variances PASS',out['variance_domain']['maximum_upper']['display_only'],flush=True)
 if not args.arithmetic_only:
  out['old_probability_gate']=old_probability_gate();print('Old LP at new P0 PASS',flush=True)
  out['endpoint']=endpoint();print('Endpoint PASS',flush=True)
 out['status']='ARITHMETIC_ONLY_PASS' if args.arithmetic_only else 'SELECTED_MOMENT_CERTIFICATE_COMPLETE_PENDING_ROOT_REVIEW'
 deps=[Path(__file__),HERE/'SELECTED_MOMENT_TRACE_PROOF.md',PRIOR/'SIZE_TRACE_CERTIFICATE.json',
       PRIOR/'certify_trace_concentration.py',OLD.HELPER,ROOT/'agents/analytic_mass_round3/MASS_MOMENT_CERTIFICATE.json',
       ROOT/'agents/joint_analytic/HALF_MODULUS_ENDPOINT_CERTIFICATE.json']
 out['dependencies']=[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in deps]
 out['external_theorem']='Independent nonidentical Berry-Esseen with C=0.56; no numerical recreation of its original proof.'
 out['seconds']=time.monotonic()-start
 out['scope']='Arb enclosures plus exact endpoint gates, conditional on reviewed general proof and inherited geometry/moments. Not all-prime connectivity or full formal verification.'
 name='SELECTED_MOMENT_ARITHMETIC.json' if args.arithmetic_only else 'SELECTED_MOMENT_CERTIFICATE.json'
 (HERE/name).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'seconds':out['seconds'],'certificate':name},indent=2),flush=True)

if __name__=='__main__':main()
