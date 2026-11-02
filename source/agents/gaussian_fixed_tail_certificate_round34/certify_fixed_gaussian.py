#!/usr/bin/env python3
"""Fixed all-frequency/all-variance certificate, Arb 192 bits.

Algorithmic adaptation of round28, with all helpers copied/reviewed here;
no predecessor module is imported and no predecessor output is modified.
No floating-point quantity is used in a mathematical comparison.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, hashlib, json, math, sys, time, traceback

sys.path.insert(0, '/tmp/baragar_math_deps')
import flint
from flint import arb, acb, ctx
ctx.prec = 192
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ETA=F(1,18); BETA=F(80,18); EPS=F(63,1000); DELTA=F(829,9180)
LAM=F(43,1000); NU=F(61,10000); D=F(1693,200); RHO=F(71,2)
C=F(6); A=14000; P=10**153; X=F(P,2); GAMMA=F(189,500)
PRIME_CUTOFF=4099; OMEGA=512; VSTOP=1600
VFLOOR=F(207); LOG2_LOWER=F(693,1000); VBOXES=2128
PI=arb.pi(); SQRT2=arb(2).sqrt()

def aa(x):
    x=F(x)
    return arb(x.numerator)/x.denominator

def box(lo,hi):
    lo,hi=F(lo),F(hi)
    return arb(aa((lo+hi)/2),aa((hi-lo)/2).upper())

def pack(x):
    assert x.is_finite()
    return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(name,data):
    path=HERE/name; temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n'); temp.replace(path)

def parameters():
    return {k:str(v) for k,v in [('eta',ETA),('beta',BETA),('epsilon',EPS),
        ('delta',DELTA),('lambda',LAM),('nu',NU),('D',D),('rho',RHO),
        ('sqrt_rho_upper',C),('P',P),('X',X),('A',A),
        ('gamma_target',GAMMA),('Omega',OMEGA),('cubic_exponent_coefficient',F(1,3)),
        ('budget_only_variance_floor',D/LAM),('lattice_variance_floor',VFLOOR),
        ('lattice_log2_lower',LOG2_LOWER)]}

DEPENDENCIES=[
    'agents/gaussian_replacement_audit_round28/INDEPENDENT_GAUSSIAN_REPLACEMENT_PROOF.md',
    'agents/gaussian_replacement_audit_round28/CENTERED_SYMMETRIZATION_STRENGTHENING.md',
    'agents/gaussian_replacement_audit_round28/certify_gaussian_replacement.py',
    'agents/gaussian_replacement_audit_round28/certify_gaussian_replacement_strengthened.py',
    'agents/probability_mechanism_round6/SELECTED_MOMENT_TRACE_PROOF.md',
    'agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json',
    'root_logs/SELECTED_MOMENT_ROOT_REVIEW.json',
    'control/SELECTED_MOMENT_BOUND_ROOT_ACCEPTANCE.json',
    'agents/selected_variance_scale_round31/SELECTED_VARIANCE_GLOBAL_PROOF.md',
    'agents/selected_variance_scale_round31/SELECTED_VARIANCE_CERTIFICATE.json',
    'control/SELECTED_VARIANCE_SCALE_ROOT_ACCEPTANCE.json',
    'control/GAUSSIAN_STRENGTHENED_TRACE_ROOT_ACCEPTANCE.json',
]

def dependency_bindings():
    return [{'path':rel,'sha256':digest(ROOT/rel),'bytes':(ROOT/rel).stat().st_size} for rel in DEPENDENCIES]

def initial_gates():
    accepted=ROOT/'control/SELECTED_VARIANCE_SCALE_ROOT_ACCEPTANCE.json'
    assert digest(accepted)=='9a2b57c3f4a7a413cc9083589bbf67938e6ae23640aab2bc0cd68ffbeb0ea290'
    data=json.loads(accepted.read_text())
    assert data['status']=='ACCEPTED_GENERAL_SELECTED_VARIANCE_SCALE_NO_NEW_ENDPOINT'
    assert F(data['rho'])==RHO
    moment=ROOT/'agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json'
    review=json.loads((ROOT/'root_logs/SELECTED_MOMENT_ROOT_REVIEW.json').read_text())
    assert review['status']=='PASS_INDEPENDENT_SELECTED_MOMENT_FULL_ARITHMETIC_AND_ENDPOINT'
    assert review['certificate_sha256']==digest(moment)
    assert digest(moment)=='ca42290f5f622d7e71ab7f68041d2467350e45a0a91a519c8f08144daf8f9992'
    expression=arb(A).log()+aa(DELTA-EPS)*aa(X).log()-arb(43540).log()
    assert expression>aa(D)>0
    assert aa(RHO)<aa(C*C)
    # The accepted centered lattice lemma gives b_i >= log(q_i)*v_i/2.
    # Thus T >= log(2)*V/2, and D <= (lambda-nu*log(2)/2)*V.
    assert arb(2).log()>aa(LOG2_LOWER)
    coefficient=LAM-NU*LOG2_LOWER/2
    assert coefficient>0 and coefficient*VFLOOR<D
    out={'status':'PASS_FIXED_BUDGET_AND_ACCEPTED_VARIANCE_SCALE',
         'parameters':parameters(),'K_cap':43540,'D_expression':pack(expression),
         'D_margin':pack(expression-aa(D)),'dependencies':dependency_bindings(),
         'lattice_variance_gate':{'log2_lower':str(LOG2_LOWER),'log2':pack(arb(2).log()),
             'coefficient_upper':str(coefficient),'coefficient_times_207':str(coefficient*VFLOOR),
             'D_minus_coefficient_times_207':str(D-coefficient*VFLOOR),
             'stronger_exact_rational_floor':str(D/coefficient),'safe_variance_floor':str(VFLOOR)},
         'producer_sha256':digest(__file__),'flint_python_version':flint.__version__,
         'precision_bits':ctx.prec}
    save('FIXED_PARAMETER_GATES.json',out)
    return out

def primes_below(stop):
    sieve=bytearray(b'\x01')*stop; sieve[:2]=b'\x00\x00'
    for q in range(2,math.isqrt(stop-1)+1):
        if sieve[q]:
            for m in range(q*q,stop,q): sieve[m]=0
    return [q for q in range(2,stop) if sieve[q]]

def local_data():
    data=[]; rows=[]; count=0
    for q in primes_below(PRIME_CUTOFF):
        logq=arb(q).log(); r=(-aa(ETA)*logq).exp(); cq=arb(q-1)/q
        epsw=(-aa(EPS)*logq).exp(); infinite=1+cq*r/(1-r)
        apow=arb(1); rpow=arb(1); sumjr=arb(0); weights=[]; a=0
        while True:
            a+=1; apow*=epsw; rpow*=r; sumjr+=a*rpow
            if apow*infinite<1: break
            weights.append((apow,(apow*cq*logq*sumjr).upper()))
            count+=1; assert a<1000
        rows.append({'q':q,'first_omitted_exponent':a,'tail_guard':pack(apow*infinite)})
        data.append((q,logq,r,cq,weights))
    lq=arb(PRIME_CUTOFF).log()
    guard=(-aa(EPS)*lq).exp()*(1+(-aa(ETA)*lq).exp())
    assert guard<1
    assert len(data)==564 and count==1355
    return data,rows,guard,count

def frequency_grid():
    return ([(F(i,16384),F(i+1,16384)) for i in range(4096)]
          + [(F(i,1024),F(i+1,1024)) for i in range(256,8192)]
          + [(F(i,128),F(i+1,128)) for i in range(1024,65536)])

def frequency_jobs(chunk_size=256):
    jobs=[]; offset=0
    for denominator,first,last in [(16384,0,4096),(1024,256,8192),(128,1024,65536)]:
        for start in range(first,last,chunk_size):
            length=min(chunk_size,last-start)
            jobs.append((offset+start-first,denominator,start,length))
        offset+=last-first
    assert offset==76544
    return jobs

_LOCAL=None; _DENOMINATOR_INVERSE=None
def worker_init():
    global _LOCAL,_DENOMINATOR_INVERSE
    ctx.prec=192
    _LOCAL=local_data()[0]
    budget=arb(A).log()+aa(DELTA-EPS)*aa(X).log()
    _DENOMINATOR_INVERSE=(-budget).exp()

def frequency_chunk(job):
    offset,denominator,start,length=job
    midpoint0=aa(F(2*start+1,2*denominator)); radius=aa(F(1,2*denominator))
    step=aa(F(1,denominator)); products=[arb(1) for _ in range(length)]
    # Exact-grid recurrence removes repeated trigonometric calls only.
    # Arb propagation encloses every exact r*exp(i*midpoint*log(q)).
    for q,logq,r,cq,weights in _LOCAL:
        sine,cosine=(logq*midpoint0).sin_cos(); unit=acb(cosine,sine)*r
        sine,cosine=(logq*step).sin_cos(); rotation=acb(cosine,sine)
        for k in range(length):
            power=acb(1); partial=acb(1); best=arb(1)
            for weight,derivative in weights:
                power*=unit; partial+=cq*power
                candidate=(weight*partial.abs_upper()+derivative*radius).upper()
                best=best.max(candidate).upper()
            products[k]*=best
            if k+1<length: unit*=rotation
    result=[]
    for product in products:
        value=(product*_DENOMINATOR_INVERSE).upper().min(arb(1)).upper()
        assert 0<value<=1
        result.append(str(value.fmpq()))
    return offset,result

def euler_table(workers):
    started=time.monotonic(); data,tails,primeguard,count=local_data()
    budget=arb(A).log()+aa(DELTA-EPS)*aa(X).log()
    values=[None]*76544; done=0; jobs=frequency_jobs()
    journal=HERE/'FREQUENCY_PROGRESS.jsonl'
    with journal.open('w') as log, ProcessPoolExecutor(max_workers=workers,initializer=worker_init) as pool:
        futures=[pool.submit(frequency_chunk,j) for j in jobs]
        for fut in as_completed(futures):
            offset,result=fut.result(); assert values[offset] is None
            values[offset:offset+len(result)]=result; done+=len(result)
            record={'completed_boxes':done,'total_boxes':len(values),'last_offset':offset,
                    'last_length':len(result),'seconds':time.monotonic()-started}
            log.write(json.dumps(record)+'\n'); log.flush()
            if done%2048==0 or done==len(values): print('Euler progress',json.dumps(record),flush=True)
    assert all(v is not None for v in values)
    out={'status':'PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES',
         'parameters':parameters(),'precision_bits':ctx.prec,'prime_tail_from':PRIME_CUTOFF,
         'prime_tail_guard':pack(primeguard),'prime_count':len(data),
         'finite_nonzero_exponents':count,'local_tails':tails,
         'frequency_grid':{'segments':[
             {'from':'0','to':'1/4','denominator':16384,'boxes':4096},
             {'from':'1/4','to':'8','denominator':1024,'boxes':7936},
             {'from':'8','to':'512','denominator':128,'boxes':64512}],
             'total_boxes':len(values)},
         'frequency_upper_values':values,'budget':pack(budget),
         'method':'Directed midpoint finite sums plus rigorous per-candidate global derivative times half width.',
         'midpoint_evaluation':'Arb complex rotation along exact rational uniform-grid chunks of at most256 cells; enclosure propagated at every multiplication.',
         'workers':workers,'seconds':time.monotonic()-started,
         'producer_sha256':digest(__file__),
         'fixed_gates_sha256':digest(HERE/'FIXED_PARAMETER_GATES.json'),
         'scope':'Every saved rational upper value bounds A_full throughout its full closed frequency interval.'}
    save('EULER_FREQUENCY_CERTIFICATE.json',out)
    print('FREQUENCY_FROZEN_SHA256',digest(HERE/'EULER_FREQUENCY_CERTIFICATE.json'),flush=True)
    return out

def normal_parts(v,L):
    v=aa(v); L=L if isinstance(L,arb) else aa(L); sigma=v.sqrt()
    # erfc is essential: 1+erf at the large negative beta argument cancels.
    p=(aa(ETA)**2*v/2-aa(ETA)*L).exp()*((aa(ETA)*v-L)/(SQRT2*sigma)).erfc()/2
    q=(aa(BETA)**2*v/2+aa(BETA)*L).exp()*((aa(BETA)*v+L)/(SQRT2*sigma)).erfc()/2
    return p,q

def normal_peak(v):
    def derivative(L):
        p,q=normal_parts(v,L)
        return -aa(ETA)*p+aa(BETA)*q
    lo=F(-100); hi=F(100)
    assert derivative(lo)>0 and derivative(hi)<0
    for _ in range(48):
        mid=(lo+hi)/2; d=derivative(mid)
        if d>0: lo=mid
        elif d<0: hi=mid
        else: raise ArithmeticError('Unresolved Gaussian derivative sign')
    assert derivative(lo)>0 and derivative(hi)<0
    p,q=normal_parts(v,box(lo,hi))
    return (p+q).upper(),lo,hi

def one_sided_normal(v):
    a=aa(v).sqrt()/18
    def inverse_mills(w):
        z=aa(w)
        return (-z*z/2).exp()/((2*PI).sqrt()*((-z/SQRT2).erfc()/2))
    lo=F(-10); hi=F(10)
    assert inverse_mills(lo)>a and inverse_mills(hi)<a
    for _ in range(48):
        mid=(lo+hi)/2; value=inverse_mills(mid)
        if value>a: lo=mid
        elif value<a: hi=mid
        else: raise ArithmeticError('Unresolved BE normal root sign')
    assert inverse_mills(lo)>a and inverse_mills(hi)<a
    w=box(lo,hi)
    return ((-(a+w)**2/2).exp()/((2*PI).sqrt()*a)).upper(),lo,hi

def variance_grid():
    edges=[F(i,16) for i in range(16*207,4097)]+[F(i) for i in range(257,1601)]
    assert len(edges)==VBOXES+1 and all(a<b for a,b in zip(edges,edges[1:]))
    return list(zip(edges,edges[1:]))

def be_tail():
    return (18/(2*PI).sqrt()+aa(F(28,25)*LAM/NU))/arb(VSTOP).sqrt()

def variance_certificate(table):
    started=time.monotonic()
    assert table['parameters']==parameters()
    assert table['producer_sha256']==digest(__file__)
    assert table['status']=='PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES'
    grids=frequency_grid(); values=[aa(F(s)) for s in table['frequency_upper_values']]
    assert len(grids)==len(values)==76544
    low=[]; highG_near=[]; highA=arb(0); far_mass=arb(0)
    for (lo,hi),aval in zip(grids,values):
        l=aa(lo); u=aa(hi)
        mass=(aa(hi-lo)*aa(ETA+BETA)/(PI*((aa(ETA)**2+l*l)*(aa(BETA)**2+l*l)).sqrt())).upper()
        if hi<=F(1,4):
            low.append((l*l/2,u**3/3,aa(RHO)*u*u/2,u**3/6+aa(C)*u**4/8,aval,mass))
        else:
            highA+=mass*aval
            if lo<1: highG_near.append((l*l/2,mass))
            else: far_mass+=mass
    assert len(low)==4096 and len(highG_near)==768
    highA=highA.upper(); far_mass=far_mass.upper(); floor=VFLOOR
    farG=(far_mass*(-aa(floor)/2).exp()).upper()
    tail=(aa(ETA+BETA)/(OMEGA*PI)*(1+(-aa(floor)*OMEGA**2/2).exp())).upper()
    critical=3*D/LAM; maximum=arb(0); worst=None; rows=[]; usedBE=0
    with (HERE/'VARIANCE_PROGRESS.jsonl').open('w') as journal:
        for index,(vl,vu) in enumerate(variance_grid()):
            v=aa(vl); tmax=aa((LAM*vu-D)/NU)
            assert tmax>=0
            peak,rootlo,roothi=normal_peak(vl)
            total=arb(0)
            for l2,u3,rhou2,pref,aval,mass in low:
                exponent_v=-v*l2
                gaussian=exponent_v.exp().upper()
                cap=(aval+gaussian).upper()
                error=(tmax*pref*(exponent_v+tmax*u3+rhou2).exp()).upper().min(cap).upper()
                total+=mass*error
            nearG=sum((mass*(-v*l2).exp().upper() for l2,mass in highG_near),arb(0)).upper()
            finite=(total.upper()+highA+nearG+farG).upper()
            replacement=(peak+finite+tail).upper()
            normalH,wlo,whi=one_sided_normal(vl)
            argmax=min(vu,max(vl,critical))
            beerror=(aa(F(28,25)/NU)*aa(LAM*argmax-D)/(aa(argmax)**3).sqrt()).upper()
            be=(normalH+beerror).upper(); full=replacement.min(be).upper()
            chosen='BE' if be<replacement else 'Gaussian_replacement'
            if chosen=='BE': usedBE+=1
            row={'V_interval':[str(vl),str(vu)],'T_upper':pack(tmax),
                 'Gaussian_root_interval':[str(rootlo),str(roothi)],'Gaussian_peak_upper':pack(peak),
                 'finite_low_error_upper':pack(total.upper()),'finite_near_high_G_upper':pack(nearG),
                 'finite_error_upper':pack(finite),'replacement_upper':pack(replacement),
                 'BE_normal_root_interval':[str(wlo),str(whi)],'BE_normal_H_upper':pack(normalH),
                 'BE_error_argmax':str(argmax),'BE_error_upper':pack(beerror),
                 'BE_upper':pack(be),'chosen_complete_bound':chosen,'total_upper':pack(full)}
            rows.append(row)
            if full>maximum: maximum=full; worst=[str(vl),str(vu)]
            record={'completed_boxes':index+1,'total_boxes':VBOXES,'V_interval':row['V_interval'],
                    'maximum_upper':str(maximum.fmpq()),'seconds':time.monotonic()-started}
            journal.write(json.dumps(record)+'\n'); journal.flush()
            if not full<aa(GAMMA):
                save('VARIANCE_TARGET_FAILURE.json',{'status':'FAILED_FIXED_GAMMA','row':row,
                     'parameters':parameters(),'frequency_sha256':digest(HERE/'EULER_FREQUENCY_CERTIFICATE.json')})
                raise ArithmeticError('Fixed gamma target failed; no certificate emitted')
            if index%100==0: print('Variance progress',index+1,'/',VBOXES,'; max',str(maximum),'seconds',round(time.monotonic()-started,2),flush=True)
    infinity=be_tail().upper(); assert infinity<aa(GAMMA)
    globalupper=maximum.max(infinity).upper(); assert globalupper<aa(GAMMA)
    out={'status':'PASS_FIXED_GAUSSIAN_TRACE_ALL_PRIMES_EXPONENTS_FREQUENCIES_VARIANCES_PENDING_ROOT_REVIEW',
         'parameters':parameters(),'precision_bits':ctx.prec,
         'variance_floor':str(floor),'variance_boxes':len(rows),
         'variance_grid':{'fine_until':'256','fine_step':'1/16','first_interval':['207','3313/16'],
                          'coarse_until':'1600','coarse_step':'1'},
         'maximum_box':worst,'maximum_upper':pack(maximum),'global_upper':pack(globalupper),
         'gamma_margin':pack(aa(GAMMA)-globalupper),'boxes_using_BE':usedBE,
         'finite_high_frequency_A_upper':pack(highA),
         'finite_high_frequency_G_treatment':{'exact_rectangles':'[1/4,1],768 cells per variance box',
             'far_mass_upper':pack(far_mass),'far_uniform_G_upper':pack(farG),
             'far_bound':'sum rectangular Fourier masses on [1,512] times exp(-Vfloor/2)'},
         'frequency_infinite_tail_upper':pack(tail),
         'infinite_variance_tail_from':VSTOP,'infinite_variance_tail_upper':pack(infinity),
         'variance_rows':rows,'all_variances_covered':True,'all_frequencies_covered':True,
         'all_primes_and_actual_exponents_covered':True,
         'euler_table_sha256':digest(HERE/'EULER_FREQUENCY_CERTIFICATE.json'),
         'producer_sha256':digest(__file__),'fixed_gates_sha256':digest(HERE/'FIXED_PARAMETER_GATES.json'),
         'dependencies':dependency_bindings(),'new_external_theorems':0,
         'inherited_external_theorem':'Independent nonidentically distributed Berry-Esseen C=0.56, already primary-source audited and root accepted.',
         'seconds':time.monotonic()-started,
         'scope':'Probability estimate only; no geometric endpoint or change to CURRENT is asserted.'}
    save('FIXED_GAUSSIAN_TRACE_CERTIFICATE.json',out)
    print(json.dumps({k:out[k] for k in ['status','maximum_box','maximum_upper','global_upper','gamma_margin','boxes_using_BE','seconds']},indent=2),flush=True)
    print('CERTIFICATE_FROZEN_SHA256',digest(HERE/'FIXED_GAUSSIAN_TRACE_CERTIFICATE.json'),flush=True)
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['all','euler','variance'],default='all')
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=8
    initial_gates()
    if args.stage in ['all','euler']: table=euler_table(args.workers)
    else: table=json.loads((HERE/'EULER_FREQUENCY_CERTIFICATE.json').read_text())
    if args.stage in ['all','variance']: variance_certificate(table)

if __name__=='__main__':
    try: main()
    except BaseException as exc:
        save('PRODUCER_FAILURE.json',{'status':'FAILED_NO_RESULT_CREDITED','exception_type':type(exc).__name__,
             'message':str(exc),'traceback':traceback.format_exc(),'producer_sha256':digest(__file__)})
        raise
