#!/usr/bin/env python3
"""All-frequency/all-variance Arb certificate for one fixed target.

The new object is the difference between the complete actual sum and the
same sum with ONLY the selected variables replaced by one matching normal.
This program imports no prior numerical producer or certificate helpers.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, sys, time

sys.path.insert(0, '/tmp/baragar_math_deps')
from flint import arb, acb, ctx
ctx.prec = 160
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ETA = F(1,18); BETA = F(5,9); EPS = F(63,1000)
DELTA = F(829,9180); LAM = F(43,1000); NU = F(61,10000)
D = F(8716,1000); C = LAM/NU; P = 10**153; X = F(P,2)
GAMMA = F(43,100); A = 18000; PRIME_CUTOFF = 4099

def aa(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

def box(lo,hi):
    return arb(aa((F(lo)+F(hi))/2), aa((F(hi)-F(lo))/2).upper())

def pack(x):
    assert x.is_finite()
    return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}

def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(name,data):
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')

def primes_below(qmax):
    sieve=bytearray(b'\x01')*qmax
    sieve[:2]=b'\x00\x00'
    for q in range(2,math.isqrt(qmax-1)+1):
        if sieve[q]:
            for m in range(q*q,qmax,q): sieve[m]=0
    return [q for q in range(2,qmax) if sieve[q]]

def parameters():
    return {k:str(v) for k,v in [('eta',ETA),('beta',BETA),('epsilon',EPS),
        ('delta',DELTA),('lambda',LAM),('nu',NU),('D',D),('c',C),
        ('P',P),('X',X),('A',A),('gamma_target',GAMMA)]}

def frequency_grid():
    # Entire [0,8], with a finer exact grid for the replacement transition.
    return [(F(i,8192),F(i+1,8192)) for i in range(2048)] + \
           [(F(i,1024),F(i+1,1024)) for i in range(256,8192)]

def local_data():
    rows=[]; data=[]; count=0
    for q in primes_below(PRIME_CUTOFF):
        logq=arb(q).log(); r=(-aa(ETA)*logq).exp()
        cq=arb(q-1)/q; epsw=(-aa(EPS)*logq).exp()
        inf=1+cq*r/(1-r)
        apow=arb(1); rpow=arb(1); sumjr=arb(0); weights=[]; a=0
        while True:
            a+=1; apow*=epsw; rpow*=r; sumjr+=a*rpow
            if apow*inf<1: break
            derivative=(apow*cq*logq*sumjr).upper()
            weights.append((apow,derivative))
            count+=1
            assert a<1000
        assert apow*inf<1
        rows.append({'q':q,'first_omitted_exponent':a,'tail_guard':pack(apow*inf)})
        data.append((q,logq,r,cq,weights))
    tail=(-aa(EPS)*arb(PRIME_CUTOFF).log()).exp()*(1+(-aa(ETA)*arb(PRIME_CUTOFF).log()).exp())
    assert tail<1
    assert len(data)==564 and count==1355
    return data,rows,tail,count

def euler_table():
    started=time.monotonic(); data,tailrows,primetail,count=local_data()
    budget=arb(A).log()+aa(DELTA-EPS)*aa(X).log()
    grid=frequency_grid(); frequency_rows=[]
    for index,(lo,hi) in enumerate(grid):
        omega=aa((lo+hi)/2); radius=aa((hi-lo)/2); product=arb(1)
        for q,logq,r,cq,weights in data:
            angle=logq*omega
            unit=acb(angle.cos(),angle.sin())*r
            power=acb(1); partial=acb(1); best=arb(1)
            for weight,derivative in weights:
                power*=unit; partial+=cq*power
                value=weight*partial.abs_upper()+derivative*radius
                best=best.max(value.upper())
            product*=best.upper()
        bound=(product*(-budget).exp()).upper().min(arb(1)).upper()
        assert 0<bound<=1
        frequency_rows.append(str(bound.fmpq()))
        if index%1000==0:
            print('Euler frequency boxes',index,'/',len(grid),'elapsed',round(time.monotonic()-started,2),flush=True)
    out={'status':'PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES',
         'parameters':parameters(),'precision_bits':ctx.prec,
         'prime_tail_from':PRIME_CUTOFF,'prime_tail_guard':pack(primetail),
         'prime_count':len(data),'finite_nonzero_exponents':count,'local_tails':tailrows,
         'frequency_grid':{'first':[0,'1/4'],'first_denominator':8192,
                           'second':['1/4',8],'second_denominator':1024,
                           'total_boxes':len(grid)},
         'frequency_upper_values':frequency_rows,'budget':pack(budget),
         'seconds':time.monotonic()-started,
         'method':'Directed midpoint finite sums plus rigorous per-candidate global derivative times half width.',
         'scope':'Each value bounds A_full throughout the full closed frequency interval; no point sampling assertion.'}
    save('EULER_FREQUENCY_CERTIFICATE.json',out)
    return out

PI=arb.pi(); SQRT2=arb(2).sqrt()

def normal_parts(v,L):
    v=aa(v); L=aa(L) if not isinstance(L,arb) else L; sigma=v.sqrt()
    p=(aa(ETA)**2*v/2-aa(ETA)*L).exp()*((aa(ETA)*v-L)/(SQRT2*sigma)).erfc()/2
    q=(aa(BETA)**2*v/2+aa(BETA)*L).exp()*((aa(BETA)*v+L)/(SQRT2*sigma)).erfc()/2
    return p,q

def normal_peak(v):
    # A/B is strictly increasing because Phi(z)/phi(z) is increasing.
    def der(L):
        p,q=normal_parts(v,L)
        return -aa(ETA)*p+aa(BETA)*q
    lo=F(-100); hi=F(100)
    assert der(lo)>0 and der(hi)<0
    for j in range(45):
        mid=(lo+hi)/2; d=der(mid)
        if d>0: lo=mid
        elif d<0: hi=mid
        else: raise ArithmeticError('Unresolved Gaussian derivative; raise precision')
    assert der(lo)>0 and der(hi)<0
    p,q=normal_parts(v,box(lo,hi))
    return (p+q).upper(),lo,hi

def be_tail(v):
    return (18/(2*PI).sqrt()+aa(F(28,25)*LAM/NU))/aa(v).sqrt()

def variance_certificate(table):
    started=time.monotonic()
    assert table['parameters']==parameters()
    assert table['status']=='PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES'
    grids=frequency_grid(); avals=[aa(F(x)) for x in table['frequency_upper_values']]
    assert len(avals)==len(grids)
    # Weight is positive and decreasing on nonnegative frequencies.
    cells=[]
    for (lo,hi),aval in zip(grids,avals):
        weight=(aa(ETA+BETA)/(PI*((aa(ETA)**2+aa(lo)**2)*(aa(BETA)**2+aa(lo)**2)).sqrt())).upper()
        cells.append((lo,hi,aval,weight*aa(hi-lo)))
    lowcells=cells[:2048]; highcells=cells[2048:]
    highA=sum((aval*mass for _,_,aval,mass in highcells),arb(0)).upper()
    variance_floor=D/LAM
    vstart=math.floor(variance_floor); vend=1600
    tail=((aa(ETA+BETA)/(8*PI))*(1+(-aa(variance_floor)*32).exp())).upper()
    rows=[]; maxupper=arb(0); maxrow=None
    for vl in range(vstart,vend):
        v0=max(F(vl),variance_floor); v1=F(vl+1)
        if v0>=v1: continue
        tmax=aa((LAM*v1-D)/NU)
        assert tmax>=0
        normal,lroot,rroot=normal_peak(v0)
        total=arb(0)
        for lo,hi,aval,mass in lowcells:
            l=aa(lo); u=aa(hi)
            gaussian=(-aa(v0)*l*l/2).exp().upper()
            cap=(aval+gaussian).upper()
            prefactor=tmax*(u**3/6+aa(C)*u**4/8)
            exponent=-aa(v0)*l*l/2+2*tmax*u**3/3+aa(C)**2*u*u/2
            # This bound applies on the entire (V,omega) rectangle.
            error=(prefactor*exponent.exp()).upper().min(cap).upper()
            total+=mass*error
        highG=sum((mass*(-aa(v0)*aa(lo)**2/2).exp().upper() for lo,_,_,mass in highcells),arb(0)).upper()
        errorfinite=total.upper()+highA+highG
        full=normal+errorfinite+tail
        if not full<aa(GAMMA):
            failure={'status':'FAILED_TARGET_GAMMA','V_interval':[str(v0),str(v1)],
                     'normal':pack(normal),'finite_error':pack(errorfinite),
                     'tail':pack(tail),'total':pack(full),'highA':pack(highA)}
            save('VARIANCE_TARGET_FAILURE.json',failure)
            raise ArithmeticError(json.dumps(failure))
        if full.upper()>maxupper:
            maxupper=full.upper(); maxrow=[str(v0),str(v1)]
        rows.append({'V_interval':[str(v0),str(v1)],'T_upper':pack(tmax),
            'normal_root_interval':[str(lroot),str(rroot)],'normal_peak_upper':pack(normal),
            'finite_error_upper':pack(errorfinite),'total_upper':pack(full.upper())})
        if vl%100==0:
            print('Variance boxes through',vl,'max',str(maxupper),'elapsed',round(time.monotonic()-started,2),flush=True)
    infinitytail=be_tail(vend)
    assert infinitytail<aa(GAMMA)
    out={'status':'PASS_NEW_GAUSSIAN_REPLACEMENT_UNIFORM_TRACE_PENDING_ROOT_REVIEW',
         'parameters':parameters(),'precision_bits':ctx.prec,'variance_floor':str(variance_floor),
         'frequency_infinite_tail_upper':pack(tail),'high_frequency_finite_A_integral_upper':pack(highA),
         'variance_box_count':len(rows),'variance_box_width':1,
         'maximum_box':maxrow,'maximum_upper':pack(maxupper),'variance_rows':rows,
         'infinite_variance_tail_from':vend,'infinite_variance_tail_upper':pack(infinitytail),
         'infinite_variance_tail_theorem':'Previously accepted nonidentical Berry-Esseen C=0.56; H<=18/sqrt(2*pi*V)',
         'all_variances_covered':True,'all_frequencies_covered':True,
         'euler_table_sha256':digest(HERE/'EULER_FREQUENCY_CERTIFICATE.json'),
         'producer_sha256':digest(__file__),'seconds':time.monotonic()-started,
         'scope':'Analytic probability estimate only. Geometry, endpoint combination and acceptance belong to root.'}
    save('GAUSSIAN_REPLACEMENT_CERTIFICATE.json',out)
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['all','euler','variance'],default='all')
    args=parser.parse_args()
    if args.stage in ['all','euler']:
        table=euler_table()
    else:
        table=json.loads((HERE/'EULER_FREQUENCY_CERTIFICATE.json').read_text())
    if args.stage in ['all','variance']:
        result=variance_certificate(table)
        print(json.dumps({k:result[k] for k in ['status','maximum_box','maximum_upper','seconds']},indent=2),flush=True)

if __name__=='__main__': main()
