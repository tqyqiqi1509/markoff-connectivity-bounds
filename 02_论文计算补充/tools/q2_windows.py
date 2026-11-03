#!/usr/bin/env python3
"""Exact q2 count/mass windows and the complete exponent tail."""
from pathlib import Path
from fractions import Fraction as F
import argparse,bisect,gzip,hashlib,json,resource,sys,time
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
resource.setrlimit(resource.RLIMIT_CPU,(90,90))
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
FILES={k:ROOT/'objects'/v for k,v in {'parameters': '8301b5dd4246bdc9bfe654f22c8130cbc3eea54fe0cfa5c3e0dc0a99189bedc3', 'positive': '1e5907eb2f407a48b40b4b824baa27fd8a35a4c7aab10519cf4e44603cbfe413', 'mass': 'fba56667a0f2ba6f0607223505afad4fbd158501f6d27178c3b9a4a38c4d5452', 'node_metadata': 'c00573d85a4b1423e0ce93d3afef60c8db2c94a27ccfff20a31b7fe0f351e32e', 'prefixes': 'f59af1c53c23a58b6a40fb1725415cf2e570dedc15cc069c109988dc96a697a3'}.items()}
BITS=160;S=1<<BITS;PAIR=1<<25;HS=1<<80
cd=lambda a,b:-(-a//b)
rat=lambda v:F(int(v['numerator']),int(v['denominator']))
up=lambda v:cd(v.numerator*S,v.denominator)
mul=lambda a,b:cd(a*b,S)
pack=lambda v:dict(numerator=str(v.numerator),denominator=str(v.denominator))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def iroot(n,k):
 if n<2:return n
 x=1<<cd(n.bit_length(),k)
 while True:
  y=((k-1)*x+n//x**(k-1))//k
  if y>=x:break
  x=y
 while (x+1)**k<=n:x+=1
 while x**k>n:x-=1
 assert x**k<=n<(x+1)**k
 return x

def bpower(B):
 t=B*B*S**13;x=iroot(t,15);y=x+(x**15!=t)
 assert y**15>=t and (y==0 or (y-1)**15<t)
 return y

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 started=time.monotonic();h=json.loads(FILES['parameters'].read_text());p=json.loads(FILES['positive'].read_text());mass=json.loads(FILES['mass'].read_text());nm=json.loads(FILES['node_metadata'].read_text())
 assert p['parameters_sha256']==sha(FILES['parameters'])
 assert nm['prefix_gzip_sha256']==sha(FILES['prefixes'])
 raw=gzip.decompress(FILES['prefixes'].read_bytes());assert hashlib.sha256(raw).hexdigest()==nm['prefix_raw_sha256'];pre=json.loads(raw)
 par=h['parameters'];meta=h['setup'];pr={v['q']:v for v in p['recursion_rows']};mr={v['q']:v for v in mass['recursion_rows']};rows={v['q']:v for v in par['active_primes']}
 assert (par['sigma'],par['epsilon_b'],par['eta'],par['b'])==('1/17','829/8100','1/18','2/15')
 assert mass['fixed_bits']==BITS
 assert meta['column_first_primes'][1]==5
 panels=[v['panels'] for v in pre];assert panels==meta['needed_prefix_panels'];assert panels[-1]==500000 and meta['density']==10000
 roots=[]
 for q in [2,3]:
  row=rows[q]
  for key,k in [('r',17),('s',18)]:
   lo,hi=rat(row[key+'_lower']),rat(row[key+'_upper']);assert lo**k<=F(1,q)<=hi**k
   roots.append(dict(q=q,kind=key,degree=k,passed=True))
  wu=rat(row['w_upper']);assert 0<wu<1 and wu**8100*F(q)**829>=1
  roots.append(dict(q=q,kind='w_upper',degree=8100,passed=True))
 piL,piU=rat(par['pi_lower']),rat(par['pi_upper']);assert 3<piL<piU<F(22,7)
 K5=up(rat(meta['Euler_tail_caps'][1]));C5m=up(rat(pr[5]['C_upper']));C5e=int(mr[5]['C_mass_scaled']);Q=rat(meta['curvature'][1]['Q_mellin_upper'])
 windows=[]
 def J(a,kind):
  if not a:return 0
  u=F(1,17) if kind=='m' else F(18,17);row=rows[2]
  lo,hi=rat(row['r_lower'])**a,rat(row['r_upper'])**a
  if kind=='E':lo/=2**a;hi/=2**a
  AU=(1-lo)/u;BU=1+hi;cc=piU*u/300;assert cc<F(1,2)
  phi=min(F(1),(1+hi)/2+(1-hi)*cc)
  jj=bisect.bisect_left(panels,BU/AU*10000);indices={max(1,min(len(pre)-1,jj)),max(1,min(len(pre)-1,jj-1))}
  vals=[]
  for j in sorted(indices):
   low=F(pre[j]['pair_sums'][1],10000*PAIR)
   harm=F(pre[-1]['harmonic_upper_scaled'][1]-pre[j]['harmonic_lower_scaled'][1],HS*PAIR)
   total=F(pre[-1]['pair_sums'][1],10000*PAIR)
   value=Q*(AU*low/piL+BU*harm/piL+2*phi*total/50)
   vals.append((value,j))
  value,j=min(vals);iv=min(K5,up(value))
  windows.append(dict(a=a,kind=kind,split_panels=panels[j],upper_scaled=str(iv)))
  return iv
 row=rows[2];wu,ru,su=[up(rat(row[k+'_upper'])) for k in ['w','r','s']]
 Binf=up(1+F(1,2)*rat(row['s_upper'])/(1-rat(row['s_upper'])));Fmax=bpower(Binf)
 a_cut=0;wpa=S
 while True:
  tm=mul(mul(wpa,Fmax),K5+C5m);te=mul(mul(wpa,Fmax),K5+C5e)
  if a_cut>0 and tm<=C5m and te<=C5e:break
  a_cut+=1;wpa=mul(wpa,wu);assert a_cut<200
 J1E=J(1,'E');cm,ce=C5m,C5e;am=ae=0;co_rows=[];wp=rp=sp=B=S;geo=0
 for a in range(1,a_cut):
  geo+=rp;wp,rp,sp=mul(wp,wu),mul(rp,ru),mul(sp,su);B+=cd(sp,2);co=mul(wp,bpower(B));jm=J(a,'m');je=min(mul(geo,J1E),jm,K5)
  vm=mul(co,jm+mul(rp,C5m));ve=mul(co,je+mul(rp,C5e))
  if vm>cm:cm,am=vm,a
  if ve>ce:ce,ae=ve,a
  co_rows.append(dict(a=a,count_scaled=str(vm),mass_scaled=str(ve),coefficient_scaled=str(co),J_count_scaled=str(jm),J_mass_scaled=str(je)))
 assert cm<=10100*S and ce<=8487*S
 print(json.dumps(dict(no3_count_display=cm/S,no3_mass_display=ce/S,no3_first_omitted_a=a_cut)),flush=True)
 out=dict(status='PASS_Q2_COUNT_MASS_WINDOWS',fixed_bits=BITS,no3=dict(C_m_scaled=str(cm),C_E_scaled=str(ce),first_omitted_a2=a_cut,omitted_count_test_scaled=str(tm),omitted_mass_test_scaled=str(te),inherited_count_next_scaled=str(C5m),inherited_mass_next_scaled=str(C5e),inherited_euler_tail_scaled=str(K5),candidates=co_rows,windows=windows),verified_parameter_root_tests=roots,dependencies=[dict(role=k,path=str(v.relative_to(ROOT)),sha256=sha(v)) for k,v in FILES.items()],source_sha256=sha(Path(__file__)),seconds=time.monotonic()-started)
 args.output.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(dict(status=out['status'],q2_tail=a_cut,output=str(args.output))))
if __name__=='__main__':main()
