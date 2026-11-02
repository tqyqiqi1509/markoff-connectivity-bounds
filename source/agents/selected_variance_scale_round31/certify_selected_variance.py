#!/usr/bin/env python3
"""All primes and all exponents via proved lattice and monotone-variance tails.

No production certificate or moment producer is imported. Arithmetic backend:
python-flint / Arb. All saved bounds are exact integers divided by 10**18.
"""
from pathlib import Path
import csv, hashlib, json, math, sys, time
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb, ctx
ctx.prec=192
Q=2**21
SCALE=10**18
LAMBDA=arb(43)/1000
NU=arb(61)/10000
C=arb(430)/61
C2=arb(184900)/3721
RHO=36

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def scaled_box(x):
    y=x*SCALE
    return int(y.lower().fmpq().floor()),int(y.upper().fmpq().ceil())
def exact_box(x):
    return {'lower':str(x.lower().fmpq()),'upper':str(x.upper().fmpq()),'display_only':str(x)}
def primes_below(n):
    sieve=bytearray(b'\1')*n
    sieve[0:2]=b'\0\0'
    for p in range(2,math.isqrt(n-1)+1):
        if sieve[p]:sieve[p*p:n:p]=b'\0'*len(range(p*p,n,p))
    return [p for p in range(2,n) if sieve[p]]

start=time.monotonic()
assert arb(Q).log()>2*C
prime_list=primes_below(Q)
prime_digest=hashlib.sha256(','.join(map(str,prime_list)).encode()).hexdigest()
moment_path=HERE/'ALL_FINITE_MOMENTS.csv'
stop_path=HERE/'ALL_PRIME_STOPS.csv'
hist={}
selected_count=excluded_count=uncertain_count=moment_count=0
max_record=None
max_scaled=-1
minimum_excluded_margin=None
minimum_selected_margin=None
max_stop_a=0
with moment_path.open('w',newline='') as mf,stop_path.open('w',newline='') as sf:
    mw=csv.writer(mf);sw=csv.writer(sf)
    mw.writerow(['q','a','v_lower_scaled','v_upper_scaled','b_lower_scaled','b_upper_scaled','W_lower_scaled','W_upper_scaled','selection'])
    sw.writerow(['q','stop_a','stop_v_lower_scaled','stop_v_upper_scaled','selected_a','uncertain_a','max_possible_selected_v_upper_scaled'])
    for pi,q in enumerate(prime_list):
        h=arb(q).log();h2=h*h;h3=h2*h
        r=(-h/18).exp();cq=arb(q-1)/q
        weights=[arb(1)]
        total=arb(1);s1=arb(0);s2=arb(0);rpower=arb(1)
        a=0;selected=[];uncertain=[];local_max=0
        while True:
            a+=1
            rpower*=r
            w=cq*rpower
            weights.append(w)
            total+=w;s1+=a*w;s2+=a*a*w
            mu=s1/total
            v=(s2/total-mu*mu)*h2
            vlo,vhi=scaled_box(v)
            # A verified stop removes this and every later exponent.
            if v>=C2:
                assert vlo*3721>=184900*SCALE
                sw.writerow([q,a,vlo,vhi,'|'.join(map(str,selected)),'|'.join(map(str,uncertain)),local_max])
                hist[a]=hist.get(a,0)+1
                max_stop_a=max(max_stop_a,a)
                break
            b=sum((wj*abs(arb(j)-mu)**3 for j,wj in enumerate(weights)),arb(0))/total*h3
            W=LAMBDA*v-NU*b
            blo,bhi=scaled_box(b);wlo,whi=scaled_box(W)
            if W<=0:
                kind='excluded';excluded_count+=1
                assert whi<=0
                margin=-whi
                if minimum_excluded_margin is None or margin<minimum_excluded_margin['margin_scaled']:
                    minimum_excluded_margin={'q':q,'a':a,'margin_scaled':margin}
            else:
                if W>0 and wlo>0:
                    kind='selected';selected.append(a);selected_count+=1
                    if minimum_selected_margin is None or wlo<minimum_selected_margin['margin_scaled']:
                        minimum_selected_margin={'q':q,'a':a,'margin_scaled':wlo}
                else:
                    kind='uncertain_treated_as_selected';uncertain.append(a);uncertain_count+=1
                assert v<RHO and vhi<RHO*SCALE
                local_max=max(local_max,vhi)
                if vhi>max_scaled:
                    max_scaled=vhi
                    max_record={'q':q,'a':a,'selection':kind,'variance':exact_box(v),'third_absolute_centered':exact_box(b),'selection_penalty':exact_box(W),'v_upper_scaled':vhi}
            mw.writerow([q,a,vlo,vhi,blo,bhi,wlo,whi,kind])
            moment_count+=1
            # No arbitrary cutoff is used. This assertion only stops accidental bugs.
            assert a<100000
        if (pi+1)%20000==0:
            print(json.dumps({'primes_done':pi+1,'q':q,'moments':moment_count,'max_q_a':None if max_record is None else [max_record['q'],max_record['a']],'max_v_upper_scaled':max_scaled,'elapsed_seconds':time.monotonic()-start}),flush=True)

assert moment_count==selected_count+excluded_count+uncertain_count
summary={
    'status':'PASS_ALL_PRIMES_ALL_EXPONENTS_SELECTED_VARIANCE_LT_36',
    'created_utc':datetime.now(timezone.utc).isoformat(),
    'parameters':{'eta':'1/18','lambda':'43/1000','nu':'61/10000','c':'430/61','c_squared':'184900/3721','rho':'36'},
    'arithmetic':{'backend':'python-flint / Arb','precision_bits':ctx.prec,'saved_fixed_denominator':str(SCALE),'selection_uncertainty_rule':'Treat every unproved exclusion as possibly selected'},
    'large_prime_tail':{'Q':Q,'log_Q':exact_box(arb(Q).log()),'twice_c':exact_box(2*C),'strict_log_Q_gt_twice_c':True,'proof':'centered lattice b>=log(q)*v/2'},
    'finite_prime_coverage':{'first_prime':prime_list[0],'last_prime':prime_list[-1],'prime_count':len(prime_list),'prime_sequence_comma_separated_sha256':prime_digest,'all_primes_q_lt_Q':True},
    'exponent_coverage':{'a_zero':'v=b=0, strictly unselected','method':'for each prime stop at first rigorously v>=c_squared; all later v strictly larger by decreasing-weight variance lemma','stop_histogram':dict(sorted(hist.items())),'maximum_stop_exponent':max_stop_a,'finite_moments_below_stop':moment_count,'stop_rows':len(prime_list)},
    'selection_counts':{'certified_selected':selected_count,'certified_excluded':excluded_count,'uncertain_safely_included':uncertain_count},
    'worst_possible_selected':max_record,
    'clean_rational_conclusions':{'v_less_than':'36','sum_vi_squared_le':'6*T','deleted_term_exponent_compensation_le':'18*omega^2','also_certified_v_upper':str(max_scaled)+'/'+str(SCALE)},
    'closest_excluded_penalty_to_zero':minimum_excluded_margin,
    'closest_selected_penalty_to_zero':minimum_selected_margin,
    'artifacts':[{ 'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in [moment_path,stop_path,Path(__file__)]],
    'elapsed_seconds':time.monotonic()-start,
    'limits':'No new p endpoint or frequency certificate is claimed; root independent review required. Shared Arb backend is not independent software verification.'
}
(HERE/'SELECTED_VARIANCE_CERTIFICATE.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)
