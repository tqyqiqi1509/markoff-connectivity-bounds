#!/usr/bin/env python3
"""Independent standard-library integer replay of positive H+ certificate.
No imports from the producer and no Arb dependency. Root/log/pi, all local
tails, curvature, every node byte/pair prefix/harmonic prefix and the whole
recursion are replayed. Pointwise transcendental node validity is separately
checked by check_nodes_chebyshev.py; this checker treats that stream as input.
"""
import argparse,gzip,hashlib,json,struct,sys,time
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
sys.set_int_max_str_digits(0)
S=1<<160
def cd(a,b):
 assert a>=0 and b>0
 return (a+b-1)//b
def mul(a,b):return cd(a*b,S)
def f(x):return F(int(x['numerator']),int(x['denominator']))
def u(x):return cd(x.numerator*S,x.denominator)
def l(x):return x.numerator*S//x.denominator
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as z:
  while b:=z.read(1048576):h.update(b)
 return h.hexdigest()
def rootfloor(n,k):
 x=1<<cd(n.bit_length(),k)
 while True:
  y=((k-1)*x+n//x**(k-1))//k
  if y>=x:break
  x=y
 while x**k>n:x-=1
 while (x+1)**k<=n:x+=1
 assert x**k<=n<(x+1)**k
 return x
def positive_factor(B):return rootfloor(B*B*S**13,15)+1
def atanh(z):
 zz=z*z;t=z;v=F(0)
 for i in range(140):v+=2*t/(2*i+1);t*=zz
 return v,v+2*t/((281)*(1-zz))
def atan(z,N):
 t=z;v=F(0)
 for i in range(N):v+=t/(2*i+1);t*=-(z*z)
 d=t/(2*N+1)
 return min(v,v+d),max(v,v+d)
def check_parameters(p):
 assert [p[k] for k in ['epsilon_b','epsilon_m','sigma','eta','b']]==['829/8100','9/100','1/17','1/18','2/15']
 rows=p['active_primes'];ps=[q for q in range(2,379) if all(q%d for d in range(2,isqrt(q)+1))]
 assert [r['q'] for r in rows]==ps and len(rows)==74
 a,b=atan(F(1,5),110);c,d=atan(F(1,239),32)
 assert 3<f(p['pi_lower'])<=16*a-4*d<=16*b-4*c<=f(p['pi_upper'])<4
 lo2,hi2=atanh(F(1,3));checks=0
 for row in rows:
  assert f(row['K_upper'])>=1 and f(row['L_upper'])>=0 and f(row['D2_upper'])>=0
  q=row['q']
  for key,num,den in [('w',829,8100),('r',1,17),('s',1,18)]:
   lo,hi=f(row[key+'_lower']),f(row[key+'_upper'])
   assert 0<lo<=hi<1
   assert lo.numerator**den*q**num<=lo.denominator**den
   assert hi.numerator**den*q**num>=hi.denominator**den
  power=q.bit_length()-1;low,high=atanh(F(q-2**power,q+2**power))
  assert f(row['logq_lower'])<=power*lo2+low<=power*hi2+high<=f(row['logq_upper'])
  w,r,s=[f(row[k+'_upper']) for k in ['w','r','s']]
  Binf=1+F(q-1,q)*s/(1-s);fm=f(row['Fmax_upper'])
  assert f(row['B_inf_upper'])==Binf and fm**15>=Binf*Binf
  E=row['first_discarded_exponent'];assert w**E*fm/(1-r)==f(row['discarded_tail_upper'])<1
  assert [x['e'] for x in row['local_candidates']]==list(range(1,E))
  B=F(1);rp=sp=F(1);s0=F(1);s1=s2=F(0)
  for c in row['local_candidates']:
   e=c['e'];rp*=r;sp*=s;B+=F(q-1,q)*sp;s0+=rp;s1+=e*rp;s2+=e*e*rp
   fu,cu=f(c['factor_upper']),f(c['coefficient_upper'])
   assert f(c['B_upper'])==B and fu**15>=B*B and cu>=w**e*fu
   assert cu*s0<=f(row['K_upper']) and cu*s1<=f(row['L_upper']) and cu*s2<=f(row['D2_upper'])
   checks+=1
 o=p['first_omitted_prime'];q=o['q'];assert q==379 and all(q%d for d in range(2,isqrt(q)+1))
 for key,num,den in [('w',829,8100),('r',1,17),('s',1,18)]:
  hi=f(o[key+'_upper']);assert 0<hi<1;assert hi.numerator**den*q**num>=hi.denominator**den
 cutoff=f(o['cutoff_factor_upper'])
 assert cutoff**15>=(1+f(o['s_upper']))**2
 assert f(o['test_upper'])==f(o['w_upper'])*(1+f(o['r_upper']))*cutoff<1
 return {'prime_count':len(rows),'local_candidates':checks,'all_roots_by_integer_powers':True,
  'logs_by_positive_atanh_series':True,'pi_by_Machin_series':True,'all_infinite_prime_and_exponent_tails':True}

def setup_replay(p,m):
 rows=p['active_primes'];n=len(rows);K=[u(f(x['K_upper'])) for x in rows]
 tails=[S]*(n+1)
 for j in range(n-1,-1,-1):tails[j]=mul(K[j],tails[j+1])
 assert m['columns']==n and m['column_prime_indices']==list(range(1,n+1))
 assert m['column_first_primes']==[r['q'] for r in rows[1:]]+[None]
 pi=f(p['pi_upper']);density=m['density'];rankN=m['rank_panels'];qms=[];qrs=[]
 for j in range(1,n+1):
  ll=sum((f(rows[k]['L_upper'])*f(rows[k]['logq_upper']) for k in range(j,n)),F(0))
  dd=sum((f(rows[k]['D2_upper'])*f(rows[k]['logq_upper'])**2 for k in range(j,n)),F(0))
  lr=sum((f(rows[k]['L_upper']) for k in range(j,n)),F(0))
  dr=sum((f(rows[k]['D2_upper']) for k in range(j,n)),F(0))
  am=dd+ll*ll;ar=dr+lr*lr
  c=m['curvature'][j-1]
  assert c['prime_index']==j and am<=f(c['A_mellin_upper']) and ar<=f(c['A_rank_upper'])
  assert am<8*density*density and ar*pi*pi<8*rankN*rankN
  qm=1/(1-am/F(8*density*density));qr=1/(1-ar*pi*pi/F(8*rankN*rankN))
  assert qm<=f(c['Q_mellin_upper']) and qr<=f(c['Q_rank_upper'])
  qms.append(u(qm));qrs.append(u(qr))
 thresholds={(t['q'],t['a']):t for t in m['thresholds']}
 assert len(thresholds)==len(m['thresholds'])
 univ={t['q']:t for t in m['universal_exponent_bounds']};needed={0}
 for lam in m['lambdas']:
  N=F(lam)*density;assert N.denominator==1;needed.add(int(N))
 count=0;maxN=max(needed)
 for j,row in enumerate(rows):
  q=row['q'];v=univ[q];A=v['first_omitted_exponent'];cap=f(v['Euler_tail_cap'])
  product=F(1)
  for k in range(j+1,n):product*=f(rows[k]['K_upper'])
  assert product<=cap and f(row['w_upper'])**A*f(row['Fmax_upper'])*(cap+1)==f(v['omitted_test_upper'])<=1
  rlo=f(row['r_lower']);rhi=f(row['r_upper']);rhoL=rhoU=S
  for a in range(1,A):
   rhoL=rhoL*rlo.numerator//rlo.denominator
   rhoU=cd(rhoU*rhi.numerator,rhi.denominator)
   t=thresholds[q,a];AU=17*F(S-rhoL,S);BU=F(S+rhoU,S)
   assert f(t['rho_lower'])==F(rhoL,S) and f(t['rho_upper'])==F(rhoU,S)
   assert f(t['A_upper'])==AU and f(t['B_upper'])==BU
   z=BU*density/AU;k=z.numerator//z.denominator+1;assert k==t['split_panels']
   needed.add(min(k,maxN));count+=1
 assert count==len(thresholds) and sorted(needed)==m['needed_prefix_panels']
 return K,tails,qms,qrs,thresholds,univ

def stream_replay(folder,mode,m):
 info=json.loads((folder/(mode+'.json')).read_text())
 packed=folder/info['prefix_file'];assert sha(packed)==info['prefix_gzip_sha256']
 raw=gzip.decompress(packed.read_bytes());assert hashlib.sha256(raw).hexdigest()==info['prefix_raw_sha256']
 records=json.loads(raw);pref={r['panels']:r for r in records}
 assert len(pref)==len(records)
 assert info['mode']==mode and info['node_scale_bits']==24 and info['pair_scale_bits']==25 and info['harmonic_scale_bits']==80
 assert info['density']==(m['density'] if mode=='mellin' else None)
 N=info['panels'];n=m['columns'];scale=1<<info['node_scale_bits'];HS=1<<info['harmonic_scale_bits']
 assert N==(m['rank_panels'] if mode=='rank' else max(F(x)*m['density'] for x in m['lambdas']))
 assert info['endpoint_nodes']==N+1 and info['columns']==n and info['pair_scale_bits']==info['node_scale_bits']+1
 assert sorted(pref)==([0,N] if mode=='rank' else m['needed_prefix_panels'])
 path=folder/info['node_file'];assert sha(path)==info['node_gzip_sha256'];codec=struct.Struct('<'+'Q'*n)
 sums=[0]*n;hl=[0]*n;hu=[0]*n;block=[0]*n;h=hashlib.sha256();bi=0;prev=None;first=None;prefix_checked=0
 with gzip.open(path,'rb') as stream:
  for node in range(N+1):
   raw=stream.read(codec.size);assert len(raw)==codec.size;h.update(raw);vals=codec.unpack(raw)
   assert vals[-1]==scale and all(scale<=x<2**64 for x in vals)
   assert all(vals[j]>=vals[j+1] for j in range(n-1))
   if first is None:first=vals
   for j,value in enumerate(vals):
    sums[j]+=value;block[j]+=value
    if mode=='mellin' and node>=2:
     v=prev[j]+value;hl[j]+=v*HS//(node-1);hu[j]+=cd(v*HS,node-1)
   if node in pref:
    prefix_checked+=1
    expected=[2*sums[j]-first[j]-vals[j] for j in range(n)]
    assert pref[node]['pair_sums']==expected
    assert pref[node]['harmonic_lower_scaled']==hl and pref[node]['harmonic_upper_scaled']==hu
   if node==info['blocks'][bi]['last_node']:
    assert block==info['blocks'][bi]['endpoint_sums'];block=[0]*n;bi+=1
   prev=vals
  assert stream.read(1)==b''
 assert bi==len(info['blocks']) and prefix_checked==len(pref) and sums==info['endpoint_totals']
 assert h.hexdigest()==info['node_raw_sha256'] and info['raw_bytes']==codec.size*(N+1)
 return info,pref,{'nodes':N+1,'values':(N+1)*n,'raw_sha256':h.hexdigest(),'all_bytes_crc_prefixes_blocks':True}

def recurse(p,m,K,tails,qm,qr,thresholds,univ,ri,rp,mi,mp):
 rows=p['active_primes'];n=len(rows);density=m['density'];pairscale=1<<mi['pair_scale_bits'];HS=1<<mi['harmonic_scale_bits']
 piL=l(f(p['pi_lower']));piU=f(p['pi_upper']);lams=list(map(F,m['lambdas']))
 W=[min(tails[j+1],cd(rp[ri['panels']]['pair_sums'][j]*qr[j],ri['panels']*pairscale)) for j in range(n)]
 def mellin(j,a):
  t=thresholds[rows[j]['q'],a];AU=u(f(t['A_upper']));BU=u(f(t['B_upper']));rho=f(t['rho_upper'])
  ans=None;chosen=None
  for lam in lams:
   N=int(lam*density);k=min(N,t['split_panels']);Sn=mp[N]['pair_sums'][j];Sk=mp[k]['pair_sums'][j]
   harmonic=mp[N]['harmonic_upper_scaled'][j]-mp[k]['harmonic_lower_scaled'][j] if k<N else 0
   assert harmonic>=0
   c=piU/(102*lam);phi=min(F(1),(1+rho)/2+(1-rho)*c) if c<F(1,2) else F(1)
   aterm=cd(AU*Sk*S,piL*density*pairscale)
   bterm=cd(BU*harmonic*S,piL*HS*pairscale)
   eterm=cd(2*u(phi)*Sn*lam.denominator,lam.numerator*density*pairscale)
   val=mul(qm[j],aterm+bterm+eterm)
   if ans is None or val<ans:ans,chosen=val,str(lam)
  return ans,chosen
 C=S;details=[]
 for j in range(n-1,-1,-1):
  row=rows[j];q=row['q'];w=u(f(row['w_upper']));r=u(f(row['r_upper']));s=u(f(row['s_upper']));fm=u(f(row['Fmax_upper']))
  best=C;maxa=0;chosen=None;a=0;wa=ra=sp=S;B=S
  while True:
   a+=1;wa=mul(wa,w);ra=mul(ra,r);sp=mul(sp,s);B+=cd((q-1)*sp,q)
   stop=mul(mul(wa,fm),tails[j+1]+C)
   if stop<=C:break
   assert a<univ[q]['first_omitted_exponent']
   J=tails[j+1];E=S;method=['Euler']
   for k in range(j+1,n):
    Q=rows[k]['q'];kap=1
    while Q**kap<q**a:kap+=1
    qr0=u(f(rows[k]['r_upper']));term=S;geo=S
    for z in range(1,kap):term=mul(term,qr0);geo+=term
    value=mul(E,mul(geo,W[k-1]))
    if value<J:J=value;method=['weighted Sperner',Q,kap]
    E=mul(E,K[k])
   if j+1<n:
    value,lam=mellin(j,a)
    if value<J:J=value;method=['half-line exponential mass',rows[j+1]['q'],lam]
   val=mul(mul(wa,positive_factor(B)),J+mul(ra,C))
   if val>best:best,maxa,chosen=val,a,method
  assert a<=univ[q]['first_omitted_exponent']
  details.append({'q':q,'C_next_scaled':str(C),'C_upper_scaled':str(best),'maximizing_a':maxa,
   'maximizing_method':chosen,'first_omitted_a':a,'tail_test_scaled':str(stop)})
  C=best
 return C,details,W

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--directory',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
 folder=Path(args.directory);h=json.loads((folder/'parameters.json').read_text());p=h['parameters'];m=h['setup'];start=time.monotonic()
 checks=check_parameters(p);print(json.dumps({'stage':'parameters checked','seconds':time.monotonic()-start}),flush=True)
 K,tails,qm,qr,th,un=setup_replay(p,m)
 ri,rp,rc=stream_replay(folder,'rank',m);print(json.dumps({'stage':'rank checked'}),flush=True)
 mi,mp,mc=stream_replay(folder,'mellin',m);print(json.dumps({'stage':'mellin checked'}),flush=True)
 C,details,W=recurse(p,m,K,tails,qm,qr,th,un,ri,rp,mi,mp)
 original=json.loads((folder/'positive157.json').read_text())
 assert F(C,S)<=f(original['Hplus_upper'])
 result={'status':'PASS' if C<=34300*S else 'TARGET FAILED','integer_fraction_only':True,'bits':160,
  'Hplus_scaled':str(C),'scale':str(S),'Hplus_integer_upper':cd(C,S),'target':34300,'below_producer_upper':True,
  'parameter_checks':checks,'rank_readback':rc,'mellin_readback':mc,'all_curvature_and_thresholds_recomputed':True,
  'all_rank_cuts_and_first_mellin_tail_recomputed':True,'recursion_rows':details,
  'parameter_file_sha256':sha(folder/'parameters.json'),'producer_result_sha256':sha(folder/'positive157.json'),
  'script_sha256':sha(Path(__file__)),'seconds_display_only':time.monotonic()-start}
 Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='recursion_rows'},indent=2),flush=True)
if __name__=='__main__':main()
