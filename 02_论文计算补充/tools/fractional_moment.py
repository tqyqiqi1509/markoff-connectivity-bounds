#!/usr/bin/env python3
"""Exact fractional cover moment, inheriting accepted same-a q2 bounds.
All new tests use integers/Fraction. No new spectral nodes or float decisions.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,resource,sys,time
sys.set_int_max_str_digits(0)
R=Path(__file__).resolve().parent.parent;HERE=Path(__file__).resolve().parent
resource.setrlimit(resource.RLIMIT_AS,(1024**3,)*2)
S=1<<160
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bind(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def rat(d):return F(int(d['numerator']),int(d['denominator']))
def up(x):return -(-x.numerator*S//x.denominator)
def mul(a,b):return -(-a*b//S)
def root_floor(n,k):
    lo,hi=0,1<<((n.bit_length()+k-1)//k)
    while hi-lo>1:
        m=(lo+hi)//2
        if m**k<=n:lo=m
        else:hi=m
    assert lo**k<=n<(lo+1)**k
    return lo
def beta_power(B):
    n=B*B*S**13;x=root_floor(n,15)
    return x+(x**15!=n)
def pack(q):return {'numerator':str(q.numerator),'denominator':str(q.denominator)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--q2',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.monotonic()
    jp=args.q2.resolve();old=json.loads(jp.read_text());assert old['fixed_bits']==160 and old['status']=='PASS_Q2_COUNT_MASS_WINDOWS'
    nc=old['no3'];assert nc['first_omitted_a2']==62 and [x['a'] for x in nc['candidates']]==list(range(1,62))
    records=[];Cno3=9525;Cf=26488
    for label,c,e in [('a0',int(nc['inherited_count_next_scaled']),int(nc['inherited_mass_next_scaled']))]+[(str(x['a']),int(x['count_scaled']),int(x['mass_scaled'])) for x in nc['candidates']]+[('tail>=62',int(nc['omitted_count_test_scaled']),int(nc['omitted_mass_test_scaled']))]:
        assert c*c*e<Cno3**3*S**3
        records.append({'label':label,'count_scaled':str(c),'mass_scaled':str(e),'cubic_test':True})
    hp=R/'objects/8301b5dd4246bdc9bfe654f22c8130cbc3eea54fe0cfa5c3e0dc0a99189bedc3'
    assert sha(hp)=='8301b5dd4246bdc9bfe654f22c8130cbc3eea54fe0cfa5c3e0dc0a99189bedc3'
    par=json.loads(hp.read_text())['parameters'];row=next(x for x in par['active_primes'] if x['q']==3)
    ru,su,wu=[rat(row[k+'_upper']) for k in ['r','s','w']]
    assert ru**17>=F(1,3) and su**18>=F(1,3) and wu**8100*3**829>=1 and 0<ru<1 and 0<su<1 and 0<wu<1
    cube=root_floor(3*S**3,3)+1;assert (cube-1)**3<=3*S**3<cube**3
    r,w,s=up(ru),up(wu),up(su);wp=sp=B=S;rpows=[S];q3=[]
    for exponent in range(17):
        if exponent:
            wp=mul(wp,w);sp=mul(sp,s);B+=-(-2*sp//3);rpows.append(mul(rpows[-1],r))
        b,res=divmod(exponent,2);co=mul(wp,beta_power(B))
        nonterminal=mul(cube,sum(rpows[2*k] for k in range(b)))
        terminal=rpows[2*b] if res==0 else mul(cube,rpows[2*b])
        value=Cno3*mul(co,nonterminal+terminal);assert value<Cf*S
        q3.append({'a':exponent,'b':b,'s':res,'coefficient_scaled':str(co),'sum_scaled':str(nonterminal+terminal),'fractional_bound_scaled':str(value),'passed':True,'display_only':value/S})
    f3=F(beta_power(up(1+F(2,3)*su/(1-su))),S)
    tail=Cno3*F(cube,S)*f3*wu**17/(1-ru*ru)
    assert tail<Cf
    out={'status':'PASS_EXACT_FRACTIONAL_SLICE','normalization':'A_frac=(sum cover e^(1/3))/d^(1/3), Z=d^(1/17)n^(829/8100)/B_1/18(n)^(2/15)','C_no3':Cno3,'C_fixed_cover_r3_2':Cf,'fixed_bits':160,'same_a2_holder_checks':records,'q3_finite':q3,'q3_tail_start':17,'q3_tail_upper':pack(tail),'tail_monotonic_factor':pack(wu),'q3_cube_root_upper_scaled':str(cube),'all_q2_exponents_inherited_and_jointly_checked':True,'all_q3_exponents_covered':True,'new_spectral_nodes':0,'dependencies':[bind(p) for p in [jp,hp]],'source_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start,'scope':'Moment statement for the same phase-zero upward cover only. Fractional Wronskian geometry and any endpoint remain separately pending. Old exact q2 candidate bounds and infinite-tail proof explicitly inherited.'}
    assert not args.output.exists();args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'no3':Cno3,'cover':Cf,'finite_max':max(x['display_only'] for x in q3),'tail_display':float(tail),'seconds':out['seconds']}))
if __name__=='__main__':main()
