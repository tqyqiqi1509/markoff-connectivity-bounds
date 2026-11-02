#!/usr/bin/env python3
"""Independent rigorous small-prime support partition certificate.

Reuses the already certified large-prime node tail from the 32768-node
interval-antichain audit, checks its entire serialized stream, and computes
all new small-prime subset products with Arb. The entire subsequent bound
recurrence is exact Fraction arithmetic rounded upward to a 2^-96 grid.
No root floating-point exploration module is imported.
"""
import argparse, gzip, hashlib, json, struct, sys, time
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import flint
from flint import arb, acb, ctx
sys.set_int_max_str_digits(0)
GRID=1<<96

def ceilf(x):return -((-x.numerator)//x.denominator)
def up(x):return F(ceilf(x*GRID),GRID)
def frac(d):
    m,e=int(d['mantissa']),int(d['exponent'])
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def pair(f):return {'numerator':str(f.numerator),'denominator':str(f.denominator)}
def parsepair(d):return F(int(d['numerator']),int(d['denominator']))
def aball(f):return arb(f.numerator)/f.denominator
def afrac(x):
    assert x.is_exact()
    m,e=x.man_exp();m,e=int(m),int(e)
    return F(m*(1<<e)) if e>=0 else F(m,1<<(-e))
def ceiling(x):return int(x.upper().ceil().fmpz())
def ipow(q,f):return arb(q).root(f.denominator)**(-f.numerator)
def primes(lim):return [q for q in range(2,lim+1) if all(q%a for a in range(2,isqrt(q)+1))]
def ratio_display(f):return aball(f).str(45)

def verify_locals(data):
    ctx.prec=384
    eps,sig=F(data['epsilon']),F(data['sigma'])
    rows=data['parameters']['active_primes'];omitted=data['parameters']['first_omitted_prime']
    assert [r['q'] for r in rows]==primes(omitted-1)
    for item in rows:
        q,E=item['q'],item['first_discarded_exponent']
        w,r=ipow(q,eps),ipow(q,sig)
        for key,x in [('w',w),('r',r)]:
            assert x.lower()>=aball(frac(item[key+'_lower']))
            assert x.upper()<=aball(frac(item[key+'_upper']))
        assert (w**E/(1-r)).upper()<=aball(frac(item['discarded_tail_upper']))<1
        rp,S,L=arb(1),arb(1),arb(0)
        for e in range(1,E):
            rp*=r;S+=rp;L+=e*rp
            assert (w**e*S).upper()<=aball(frac(item['K_upper']))
            assert (w**e*L).upper()<=aball(frac(item['L_upper']))
    assert (ipow(omitted,eps)*(1+ipow(omitted,sig))).upper()<1
    return {'all_local_primes_checked':len(rows),'first_inactive_prime':omitted,
            'verification_precision':384,'all_roots_K_L_and_exponent_tails_verified':True}

def integrate_subsets(data,src,out,small_count):
    info=data['integration'];rows=data['parameters']['active_primes']
    panels=info['panels'];n=len(rows);start=small_count+1
    small=rows[1:start];assert small_count>=1 and start<n
    count=1<<small_count;scale=1<<info['node_scale_bits']
    raw=src.parent/info['node_file']
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==info['node_gzip_sha256']
    oldcodec=struct.Struct('<'+'Q'*(n+1));codec=struct.Struct('<'+'Q'*count)
    assert oldcodec.size==info['row_bytes']
    dest=out.with_suffix('.subset_nodes.u64.gz')
    ctx.prec=192;eps,sig=F(data['epsilon']),F(data['sigma']);pi=arb.pi()
    prep=[]
    for item in small:
        w,r=ipow(item['q'],eps),ipow(item['q'],sig)
        prep.append((r,[w**e for e in range(item['first_discarded_exponent'])]))
    totals=[0]*count;oldtotals=[0]*(n+1);oldblock=[0]*(n+1)
    blocks=[];block=[0]*count;newhash=hashlib.sha256();oldhash=hashlib.sha256();bi=0
    started=time.monotonic()
    with gzip.open(raw,'rb') as source,gzip.open(dest,'wb',compresslevel=6) as output:
        for node in range(panels):
            chunk=source.read(oldcodec.size);assert len(chunk)==oldcodec.size
            oldhash.update(chunk);oldvalues=oldcodec.unpack(chunk)
            assert oldvalues[-1]==scale and all(a>=b for a,b in zip(oldvalues,oldvalues[1:]))
            for j,v in enumerate(oldvalues):oldtotals[j]+=v;oldblock[j]+=v
            z=acb(0,pi*(2*node+1)/(2*panels)).exp();gs=[]
            for r,wp in prep:
                term,poly,g=acb(1),acb(1),arb(1);rz=r*z
                for e in range(1,len(wp)):
                    term*=rz;poly+=term
                    g=max(g,(wp[e]*poly.abs_upper()).upper())
                gs.append(g)
            products=[arb(1)]*count;values=[oldvalues[start]]*count
            large=arb(oldvalues[start])/scale
            for mask in range(1,count):
                low=mask&-mask;bit=low.bit_length()-1
                products[mask]=(products[mask^low]*gs[bit]).upper()
                values[mask]=ceiling(products[mask]*large*scale)
                assert scale<=values[mask]<2**64
            encoded=codec.pack(*values);output.write(encoded);newhash.update(encoded)
            for j,v in enumerate(values):totals[j]+=v;block[j]+=v
            if node==info['blocks'][bi]['last_node']:
                assert oldblock==info['blocks'][bi]['tail_sums'];oldblock=[0]*(n+1);bi+=1
            if (node+1)%1024==0 or node+1==panels:
                blocks.append({'last_node':node,'subset_sums':block});block=[0]*count
                status={'stage':'subset_integrals','completed_nodes':node+1,'panels':panels,
                        'seconds_display_only':round(time.monotonic()-started,3)}
                print(json.dumps(status),flush=True);out.with_suffix('.checkpoint.json').write_text(json.dumps(status)+'\n')
        assert source.read(1)==b''
    assert oldhash.hexdigest()==info['node_raw_sha256']
    assert oldtotals==[int(t['node_sum_integer']) for t in info['tails']]
    assert bi==len(info['blocks'])
    # Full readback of every new subset value, block sum and total.
    sums=[0]*count;subblock=[0]*count;h=hashlib.sha256();bi=0
    with gzip.open(dest,'rb') as f:
        for i in range(panels):
            b=f.read(codec.size);assert len(b)==codec.size;h.update(b)
            vs=codec.unpack(b)
            for j,v in enumerate(vs):sums[j]+=v;subblock[j]+=v
            if i==blocks[bi]['last_node']:
                assert subblock==blocks[bi]['subset_sums'];subblock=[0]*count;bi+=1
        assert f.read(1)==b''
    assert sums==totals and h.digest()==newhash.digest() and bi==len(blocks)
    # Rebuild all large tail bounds with pure rational outward arithmetic.
    ctx.prec=384;piU=afrac(arb.pi().upper());D=piU/(2*panels)
    K=[frac(r['K_upper']) for r in rows];L=[frac(r['L_upper']) for r in rows]
    Ktail=[F(1)]*(n+1);correction=[F(1)]*(n+1);W=[F(1)]*(n+1)
    for j in range(n-1,-1,-1):
        Ktail[j]=up(K[j]*Ktail[j+1]);correction[j]=up((1+L[j]*D)*correction[j+1])
        W[j]=min(Ktail[j],up(F(oldtotals[j],panels*scale)*correction[j]))
    Ks=[F(1)]*count;Fs=[F(1)]*count;Ws=[F(1)]*count;subrecords=[]
    for mask in range(count):
        if mask:
            low=mask&-mask;bit=low.bit_length()-1
            Ks[mask]=up(Ks[mask^low]*K[bit+1]);Fs[mask]=up(Fs[mask^low]*(1+L[bit+1]*D))
        cap=up(Ks[mask]*Ktail[start]);corr=up(Fs[mask]*correction[start])
        Ws[mask]=min(cap,up(F(totals[mask],panels*scale)*corr))
        subrecords.append({'mask':mask,'node_sum_integer':str(totals[mask]),
                           'K_upper':pair(cap),'correction_upper':pair(corr),'W_upper':pair(Ws[mask])})
    meta={'source_file':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
          'source_node_sha256':info['node_raw_sha256'],'source_all_nodes_readback':True,
          'source_values_checked':panels*(n+1),'source_blocks_checked':len(info['blocks']),
          'small_primes':[r['q'] for r in small],'large_tail_first_prime':rows[start]['q'],
          'panels':panels,'node_scale_bits':info['node_scale_bits'],'subset_count':count,
          'node_file':dest.name,'row_bytes':codec.size,'raw_sha256':newhash.hexdigest(),
          'gzip_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
          'new_all_nodes_readback':True,'new_values_checked':panels*count,
          'pi_upper':pair(piU),'subset_integrals':subrecords,'blocks':blocks,
          'large_tail_integrals':[{'q':rows[j]['q'],'K_upper':pair(Ktail[j]),
                                  'W_upper':pair(W[j]),'correction_upper':pair(correction[j])}
                                 for j in range(start,n)],
          'seconds_display_only':time.monotonic()-started}
    return W,Ws,Ktail,Ks,meta

def support_recursion(data,W,Ws,Ktail,Ks,small_count):
    rows=data['parameters']['active_primes'];n=len(rows);start=small_count+1
    q=[r['q'] for r in rows];K=[frac(r['K_upper']) for r in rows]
    w=[frac(r['w_upper']) for r in rows];r=[frac(r['r_upper']) for r in rows]
    def step(j,C,indices,Wvalues,K0,forced=None):
        A=1
        while w[j]**A*(K0+C)>C:A+=1
        assert w[j]**A*(K0+C)<=C
        best,besta=C,0;tests=[];one=None
        for a in (range(1,A) if forced is None else [forced]):
            J,cut,kap,E=K0,None,None,F(1)
            qa=q[j]**a
            for index,Wv in zip(indices,Wvalues):
                k,Q=1,q[index]
                while Q**k<qa:k+=1
                trial=up(E*k*Wv)
                if trial<J:J,cut,kap=trial,Q,k
                E=up(E*K[index])
            value=up(w[j]**a*(J+r[j]**a*C))
            tests.append({'a':a,'J_upper':pair(J),'cut_prime':cut,'layers':kap,'value_upper':pair(value)})
            if a==1:one=value
            if value>best:best,besta=value,a
        return (one if forced is not None else best),{'q':q[j],'C_next':pair(C),'K_tail':pair(K0),
                'first_omitted_exponent':A,'tail_test_upper':pair(w[j]**A*(K0+C)),
                'best_exponent':besta,'C_upper':pair(best),'candidates':tests}
    C=F(1);big_records=[]
    for j in range(n-1,start-1,-1):
        C,rec=step(j,C,list(range(j+1,n)),W[j+1:n],Ktail[j+1]);big_records.append(rec)
    Cbig=C;count=1<<small_count;constants=[];logs=[]
    for mask in range(count):
        included=[j for j in range(1,start) if mask&(1<<(j-1))]
        C=Cbig;records=[]
        def tail(j):
            rem=[x for x in included if x>j];indices=rem+list(range(start,n));ws=[]
            for x in rem:
                sub=sum(1<<(h-1) for h in rem if h>=x);ws.append(Ws[sub])
            ws+=W[start:n]
            sub=sum(1<<(h-1) for h in rem)
            K0=up(Ks[sub]*Ktail[start])
            return indices,ws,K0
        for j in included[::-1]:
            indices,ws,K0=tail(j);C,rec=step(j,C,indices,ws,K0);records.append(rec)
        Codd=C;indices,ws,K0=tail(0)
        Call,rec=step(0,C,indices,ws,K0);records.append(rec)
        Cone,_=step(0,C,indices,ws,K0,forced=1)
        constants.append({'mask':mask,'allowed_small_primes':[q[j] for j in included],
                          'C_odd_upper':pair(Codd),'C_all_two_upper':pair(Call),'C_v2_one_upper':pair(Cone),
                          'C_all_two_display_only':ratio_display(Call),'C_v2_one_display_only':ratio_display(Cone)})
        logs.append({'mask':mask,'rows':records})
    pairs=[];full=count-1
    for mask in range(count):
        value=parsepair(constants[mask]['C_all_two_upper'])+parsepair(constants[full^mask]['C_v2_one_upper'])
        pairs.append({'high_v2_mask':mask,'low_v2_mask':full^mask,'joint_upper':pair(value),'display_only':ratio_display(value)})
    worst=max(pairs,key=lambda d:parsepair(d['joint_upper']));joint=parsepair(worst['joint_upper'])
    return ceilf(joint),{'arithmetic':'Exact Fraction with positive upward 2^-96 rounding',
         'large_tail_constant':pair(Cbig),'large_tail_rows':big_records,'support_rows':logs,
         'support_constants':constants,'all_assignments':pairs,'worst_assignment':worst,
         'joint_integer_upper':ceilf(joint),'joint_display_only':ratio_display(joint)}

def rootceil(n,k):
    lo,hi=0,1
    while hi**k<n:hi*=2
    while hi-lo>1:
        m=(lo+hi)//2
        if m**k<n:lo=m
        else:hi=m
    assert lo**k<n<=hi**k
    return hi

def endpoint(C,power):
    P=10**power;K=181265
    left=(216*(2*K+3))**900*27**420*C**6660*(P+1)**673
    right=(2*P)**900
    T=rootceil((P+1)**5,44)
    safe=36864*2541933**4*2541935**2*T**6;cage=256*2541933**4*T**4
    ratio=(arb(216)*(arb(K)+arb(3)/2)*arb(27)**(arb(7)/15)*arb(C)**(arb(37)/5)*arb(P+1)**(arb(673)/900)/P)
    return {'power':power,'C':C,'K_B':K,'T':str(T),'alpha':'673/900',
        'checks':{'exact_main':left<right,'safe8':safe<P,'cage':cage<P,'root_exact':(T-1)**44<(P+1)**5<=T**44},
        'left_sha256':hashlib.sha256(str(left).encode()).hexdigest(),'right_sha256':hashlib.sha256(str(right).encode()).hexdigest(),
        'ratio_display_only':ratio.str(45),'ratio_upper':pair(afrac(ratio.upper())),
        'safe8_left':str(safe),'cage_left':str(cage)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='baragar_run/audit/interval_antichain/arb_178_full32768.json')
    ap.add_argument('--small',type=int,default=6);ap.add_argument('--power',type=int,default=169)
    ap.add_argument('--output',default='baragar_run/audit/disjoint_support/support_169_32768.json');args=ap.parse_args()
    src,out=Path(args.input),Path(args.output);data=json.loads(src.read_text())
    assert data['epsilon']=='4/45' and data['sigma']=='1/18'
    local=verify_locals(data)
    W,Ws,Ktail,Ks,integ=integrate_subsets(data,src,out,args.small)
    C,rec=support_recursion(data,W,Ws,Ktail,Ks,args.small)
    ep=endpoint(C,args.power);low=endpoint(C,args.power-1)
    result={'status':'PASSED' if all(ep['checks'].values()) else 'ENDPOINT FAILED',
            'epsilon':'4/45','sigma':'1/18','delta':'9/100','eta':'1/18',
            'flint_version':flint.__version__,'node_Arb_precision':192,'local_verification':local,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'integration':integ,'recursion':rec,'endpoint':ep,'lower_endpoint_diagnostic':low}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'C':C,'joint':rec['joint_display_only'],
                      'worst_mask':rec['worst_assignment']['high_v2_mask'],'endpoint':ep,'lower':low},indent=2),flush=True)
    if result['status']!='PASSED':raise SystemExit(1)
if __name__=='__main__':main()
