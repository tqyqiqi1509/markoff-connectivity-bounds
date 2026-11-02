#!/usr/bin/env python3
"""Rigorous H+ producer: eps=829/8100,sigma=1/17,eta=1/18,b=2/15.
Arb only encloses local roots, finite polynomials and trig values.
All numerical decisions outside node generation use exact rational arithmetic.
All nonempty odd-prime tails plus empty tail are retained.
"""
import argparse, gzip, hashlib, json, struct, sys, time
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from flint import arb, acb, ctx
import flint
sys.set_int_max_str_digits(0)
BITS=128; SCALE=1<<BITS; ROOTBITS=192; ROOTSCALE=1<<ROOTBITS
NODEBITS=24; NS=1<<NODEBITS; HBITS=80; HS=1<<HBITS
LAMBDAS=list(map(F,['.08','.12','.18','.25','.35','.5','.75','1','1.5','2','3','4','6','8','12','18','25','35','50']))

def cd(a,b):return (a+b-1)//b
def ceilf(x):return cd(x.numerator,x.denominator)
def up(x):return F(ceilf(x*SCALE),SCALE)
def down_grid(x,s):return F(x.numerator*s//x.denominator,s)
def up_grid(x,s):return F(ceilf(x*s),s)
def pair(x):return {'numerator':str(x.numerator),'denominator':str(x.denominator)}
def rat(x):return F(int(x['numerator']),int(x['denominator']))
def afr(x):
 assert x.is_exact() and x.is_finite()
 m,e=map(int,x.man_exp())
 return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def al(x):return down_grid(afr(x.lower()),ROOTSCALE)
def au(x):return up_grid(afr(x.upper()),ROOTSCALE)
def ball(x):return arb(x.numerator)/x.denominator
def show(x):return ball(x).str(48)
def shafile(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def primes(n):return [q for q in range(2,n+1) if all(q%d for d in range(2,isqrt(q)+1))]
def positive_root_upper(x):
 # x^(2/15), with x an exact positive rational.
 u=au(ball(x).root(15)**2)
 assert u**15>=x*x
 return u
def initialize():
 ctx.prec=384
 rows=[];active=[]
 for q in primes(512):
  w=arb(q).root(8100)**(-829);r=arb(q).root(17)**(-1);s=arb(q).root(18)**(-1)
  wu,ru,su=au(w),au(r),au(s)
  cut=wu*(1+ru)*positive_root_upper(1+su)
  if cut<1:
   omitted={'q':q,'w_upper':pair(wu),'r_upper':pair(ru),'s_upper':pair(su),
     'cutoff_factor_upper':pair(positive_root_upper(1+su)),'test_upper':pair(cut)}
   break
  BinfU=1+F(q-1,q)*su/(1-su);Fmax=positive_root_upper(BinfU)
  E=1
  while wu**E*Fmax/(1-ru)>=1:E+=1
  K,L,D2=F(1),F(0),F(0)
  rp,sp,B,s0,s1,s2=F(1),F(1),F(1),F(1),F(0),F(0)
  coeff=[arb(1)];weights=[]
  spb,Bb=arb(1),arb(1)
  for e in range(1,E):
   rp*=ru;sp*=su;B+=F(q-1,q)*sp
   fu=positive_root_upper(B);cu=up(wu**e*fu)
   s0+=rp;s1+=e*rp;s2+=e*e*rp
   K=max(K,up(cu*s0));L=max(L,up(cu*s1));D2=max(D2,up(cu*s2))
   weights.append({'e':e,'B_upper':pair(B),'factor_upper':pair(fu),'coefficient_upper':pair(cu)})
   spb*=s;Bb+=(1-arb(1)/q)*spb;coeff.append(w**e*Bb.root(15)**2)
  row={'q':q,'first_discarded_exponent':E,'K_upper':pair(K),'L_upper':pair(L),
   'D2_upper':pair(D2),'B_inf_upper':pair(BinfU),'Fmax_upper':pair(Fmax),
   'discarded_tail_upper':pair(wu**E*Fmax/(1-ru)),'local_candidates':weights}
  for key,v in [('w',w),('r',r),('s',s),('logq',arb(q).log())]:
   row[key+'_lower']=pair(al(v));row[key+'_upper']=pair(au(v))
  rows.append(row);active.append({'q':q,'r':r,'coeff':coeff,'logq':arb(q).log()})
 assert len(rows)==74 and omitted['q']==379 and rows[-1]['q']==373
 p={'epsilon_b':'829/8100','epsilon_m':'9/100','sigma':'1/17','eta':'1/18','b':'2/15',
  'root_bits':ROOTBITS,'local_arithmetic_bits':BITS,'local_precision':384,
  'active_primes':rows,'first_omitted_prime':omitted,'pi_lower':pair(al(arb.pi())),'pi_upper':pair(au(arb.pi()))}
 return active,p

def setup(p,density,rankN):
 rows=p['active_primes'];n=len(rows);K=[rat(x['K_upper']) for x in rows]
 KK=[F(1)]*(n+1)
 for j in range(n-1,-1,-1):KK[j]=up(K[j]*KK[j+1])
 universal=[];thresholds=[];needed={0}
 for lam in LAMBDAS:
  N=lam*density;assert N.denominator==1;needed.add(int(N))
 for j,row in enumerate(rows):
  q=row['q'];w=rat(row['w_upper']);fm=rat(row['Fmax_upper']);A=1
  while w**A*fm*(KK[j+1]+1)>1:A+=1
  universal.append({'q':q,'first_omitted_exponent':A,'Euler_tail_cap':pair(KK[j+1]),
    'omitted_test_upper':pair(w**A*fm*(KK[j+1]+1))})
  lo=F(1);hi=F(1);rlo=rat(row['r_lower']);rhi=rat(row['r_upper'])
  for a in range(1,A):
   lo=down_grid(lo*rlo,1<<160);hi=up_grid(hi*rhi,1<<160)
   AU=17*(1-lo);BU=1+hi
   t=BU*density/AU;k=t.numerator//t.denominator+1
   needed.add(min(k,int(max(LAMBDAS)*density)))
   thresholds.append({'q':q,'a':a,'rho_lower':pair(lo),'rho_upper':pair(hi),
    'A_upper':pair(AU),'B_upper':pair(BU),'split_panels':k})
 # Node columns j=1..n, where n is the empty odd tail.
 corrections=[];curves=[];sl=sd=sr=dr=F(0);piU=rat(p['pi_upper'])
 for j in range(n,0,-1):
  if j<n:
   row=rows[j];L=rat(row['L_upper']);D=rat(row['D2_upper']);log=rat(row['logq_upper'])
   sl+=L*log;sd+=D*log*log;sr+=L;dr+=D
  Am=up(sd+sl*sl);Ar=up(dr+sr*sr)
  assert Am<8*density*density and Ar*piU*piU<8*rankN*rankN
  Qm=up(1/(1-Am/F(8*density*density)))
  Qr=up(1/(1-Ar*piU*piU/F(8*rankN*rankN)))
  curves.append({'first_prime':rows[j]['q'] if j<n else None,'prime_index':j,
   'A_mellin_upper':pair(Am),'A_rank_upper':pair(Ar),'Q_mellin_upper':pair(Qm),'Q_rank_upper':pair(Qr)})
 curves=curves[::-1]
 return {'density':density,'rank_panels':rankN,'lambdas':[str(x) for x in LAMBDAS],
  'columns':n,'column_prime_indices':list(range(1,n+1)),
  'column_first_primes':[r['q'] for r in rows[1:]]+[None],
  'column_order':'tails beginning at each odd active prime 3..373, then empty tail',
  'universal_exponent_bounds':universal,'thresholds':thresholds,
  'needed_prefix_panels':sorted(needed),'curvature':curves,'Euler_tail_caps':[pair(KK[j]) for j in range(1,n+1)]}

def produce_nodes(active,p,meta,out,mode):
 n=len(active);columns=n
 density=meta['density'];N=meta['rank_panels'] if mode=='rank' else int(max(LAMBDAS)*density)
 path=out/(mode+'.nodes.u64.gz');codec=struct.Struct('<'+'Q'*columns)
 ctx.prec=192
 needed={0,N} if mode=='rank' else set(meta['needed_prefix_panels'])
 total=[0]*columns;hl=[0]*columns;hu=[0]*columns
 endpoint_totals=[0]*columns;block=[0]*columns;blocks=[]
 records=[{'panels':0,'pair_sums':[0]*columns,'harmonic_lower_scaled':[0]*columns,'harmonic_upper_scaled':[0]*columns}]
 rawhash=hashlib.sha256();prev=None;start=time.monotonic()
 with gzip.open(path,'wb',compresslevel=3) as stream:
  for node in range(N+1):
   t=arb(node)*arb.pi()/N if mode=='rank' else arb(node)/density
   gs=[]
   for row in active[1:]:
    omega=arb(1) if mode=='rank' else row['logq']
    z=acb(0,t*omega).exp();rz=row['r']*z
    term=acb(1);poly=acb(1);g=arb(1)
    for c in row['coeff'][1:]:
     term*=rz;poly+=term;g=max(g,(c*poly.abs_upper()).upper())
    gs.append(g)
   vals=[NS]*columns;tail=arb(1)
   for j in range(len(gs)-1,-1,-1):
    tail=(tail*gs[j]).upper();vals[j]=int((tail*NS).upper().ceil().fmpz())
   assert all(NS<=v<(1<<64) for v in vals)
   assert all(vals[j]>=vals[j+1] for j in range(columns-1))
   encoded=codec.pack(*vals);stream.write(encoded);rawhash.update(encoded)
   for j,v in enumerate(vals):endpoint_totals[j]+=v;block[j]+=v
   if prev is not None:
    cell=node-1
    for j in range(columns):
     v=prev[j]+vals[j];total[j]+=v
     if mode=='mellin' and cell>0:
      numerator=v*HS;hl[j]+=numerator//cell;hu[j]+=cd(numerator,cell)
    if node in needed:
     records.append({'panels':node,'pair_sums':total.copy(),
      'harmonic_lower_scaled':hl.copy(),'harmonic_upper_scaled':hu.copy()})
   prev=vals
   if (node+1)%8192==0 or node==N:
    blocks.append({'first_node':(node//8192)*8192,'last_node':node,'endpoint_sums':block})
    block=[0]*columns
    progress={'mode':mode,'endpoint_nodes_completed':node+1,'endpoint_nodes_total':N+1,
      'seconds_display_only':round(time.monotonic()-start,3)}
    (out/(mode+'.progress.json')).write_text(json.dumps(progress)+'\n');print(json.dumps(progress),flush=True)
 assert [r['panels'] for r in records]==sorted(needed)
 pref=out/(mode+'.prefixes.json.gz')
 raw=json.dumps(records,separators=(',',':')).encode()
 with gzip.open(pref,'wb',compresslevel=6) as f:f.write(raw)
 info={'mode':mode,'endpoint_nodes':N+1,'panels':N,'density':density if mode=='mellin' else None,
  'columns':columns,'node_scale_bits':NODEBITS,'pair_scale_bits':NODEBITS+1,
  'harmonic_scale_bits':HBITS,'node_precision':192,'node_file':path.name,'node_gzip_sha256':shafile(path),
  'node_raw_sha256':rawhash.hexdigest(),'raw_bytes':codec.size*(N+1),'endpoint_totals':endpoint_totals,
  'blocks':blocks,'prefix_file':pref.name,'prefix_raw_sha256':hashlib.sha256(raw).hexdigest(),
  'prefix_gzip_sha256':shafile(pref),'prefix_count':len(records),'seconds_display_only':time.monotonic()-start}
 (out/(mode+'.json')).write_text(json.dumps(info,indent=2)+'\n')
 return info,records

def recursion(p,meta,rank_records,mellin_records):
 rows=p['active_primes'];n=len(rows);thresholds={(t['q'],t['a']):t for t in meta['thresholds']}
 pref={t['panels']:t for t in mellin_records};univ={x['q']:x for x in meta['universal_exponent_bounds']}
 density=meta['density'];pairscale=2*NS;piL=rat(p['pi_lower']);piU=rat(p['pi_upper'])
 K=[rat(x['K_upper']) for x in rows];caps=[rat(x) for x in meta['Euler_tail_caps']]
 W=[min(caps[j],up(rat(meta['curvature'][j]['Q_rank_upper'])*F(v,meta['rank_panels']*pairscale)))
   for j,v in enumerate(rank_records[-1]['pair_sums'])]
 @lru_cache(None)
 def mellin(j,a):
  tr=thresholds[rows[j]['q'],a];AU=rat(tr['A_upper']);BU=rat(tr['B_upper']);rho=rat(tr['rho_upper'])
  col=j # tail after prime index j is column j.
  best=None;choice=None
  for lam in LAMBDAS:
   N=int(lam*density);k=min(N,tr['split_panels']);end=pref[N];beg=pref[k]
   Sn=end['pair_sums'][col];Sk=beg['pair_sums'][col]
   harm=end['harmonic_upper_scaled'][col]-beg['harmonic_lower_scaled'][col] if k<N else 0
   assert harm>=0
   c=piU/(102*lam)
   phi=min(F(1),(1+rho)/2+(1-rho)*c) if c<F(1,2) else F(1)
   val=AU*F(Sk,density*pairscale)/piL+BU*F(harm,HS*pairscale)/piL+2*phi*F(Sn,density*pairscale)/lam
   val=up(rat(meta['curvature'][col]['Q_mellin_upper'])*val)
   if best is None or val<best:best,choice=val,str(lam)
  return best,choice
 details=[];C=F(1)
 for j in range(n-1,-1,-1):
  row=rows[j];q=row['q'];w=rat(row['w_upper']);r=rat(row['r_upper']);s=rat(row['s_upper'])
  Ktail=caps[j];fm=rat(row['Fmax_upper']);A=1
  while w**A*fm*(Ktail+C)>C:A+=1
  assert A<=univ[q]['first_omitted_exponent']
  stop=w**A*fm*(Ktail+C);best=C;besta=0;tests=[]
  sp=F(1);B=F(1)
  for a in range(1,A):
   sp*=s;B+=F(q-1,q)*sp;fu=positive_root_upper(B)
   J=Ktail;method=['Euler'];E=F(1)
   # All rank cuts, with an explicitly multiplied removed-prime prefix.
   for index in range(j+1,n):
    Q=rows[index]['q'];kap=1
    while Q**kap<q**a:kap+=1
    assert Q**kap>=q**a and (kap==1 or Q**(kap-1)<q**a)
    qr=rat(rows[index]['r_upper']);geo=sum((qr**h for h in range(kap)),F(0))
    val=up(E*geo*W[index-1])
    if val<J:J=val;method=['weighted Sperner',Q,kap]
    E=up(E*K[index])
   if j+1<n:
    val,lam=mellin(j,a)
    if val<J:J=val;method=['half-line exponential mass',rows[j+1]['q'],lam]
   value=up(w**a*fu*(J+r**a*C))
   tests.append({'a':a,'B_upper':pair(B),'factor_upper':pair(fu),'J_upper':pair(J),'method':method,'value_upper':pair(value)})
   if value>best:best,besta=value,a
  details.append({'q':q,'C_next_upper':pair(C),'Euler_tail_upper':pair(Ktail),'first_omitted_exponent':A,
   'omitted_tail_test_upper':pair(stop),'maximizing_exponent':besta,'C_upper':pair(best),'candidates':tests})
  C=best
 result={'status':'POSITIVE MOMENT H+ <=34300 CERTIFIED' if C<=34300 else 'TARGET FAILED',
  'Hplus_upper':pair(C),'Hplus_upper_display_only':show(C),'Hplus_integer_upper':ceilf(C),'target':34300,
  'rank_W_upper':[pair(x) for x in W],'recursion_rows':details,'mellin_cache_entries':mellin.cache_info().currsize,
  'domain':'all integers n>=1 and real d>0; all prime powers including all powers of 2'}
 return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output-dir',required=True);ap.add_argument('--density',type=int,default=10000)
 ap.add_argument('--rank-panels',type=int,default=32768);ap.add_argument('--reuse-nodes',action='store_true')
 args=ap.parse_args();out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True);active,p=initialize();meta=setup(p,args.density,args.rank_panels)
 header={'parameters':p,'setup':meta,'flint_version':flint.__version__,'source_sha256':shafile(Path(__file__)),
  'analytic_dependencies':['spectral_curvature_proof.md','HALFLINE_EXPONENTIAL_SELF_CONTAINED_PROOF.md','GENERAL_BLOCK_AND_MOMENT_PROOF.md']}
 (out/'parameters.json').write_text(json.dumps(header,indent=2)+'\n')
 print(json.dumps({'setup':'complete','active_primes':len(active),'odd_tail_columns':len(active),'needed_prefixes':len(meta['needed_prefix_panels'])}),flush=True)
 if args.reuse_nodes:
  rr=json.loads(gzip.decompress((out/'rank.prefixes.json.gz').read_bytes()))
  mr=json.loads(gzip.decompress((out/'mellin.prefixes.json.gz').read_bytes()))
 else:
  _,rr=produce_nodes(active,p,meta,out,'rank');_,mr=produce_nodes(active,p,meta,out,'mellin')
 ctx.prec=384
 result=recursion(p,meta,rr,mr)
 result['parameters_sha256']=shafile(out/'parameters.json');result['source_sha256']=shafile(Path(__file__))
 (out/'positive157.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['recursion_rows','rank_W_upper']},indent=2),flush=True)
if __name__=='__main__':main()
