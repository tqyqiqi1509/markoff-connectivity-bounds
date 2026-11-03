#!/usr/bin/env python3
"""Exact rational endpoint gates for the fractional-cover proof."""
from fractions import Fraction as Q
from pathlib import Path
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import argparse, hashlib, json, math, sys, time, traceback
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from exact_arithmetic import (fraction, file_sha, integer_digest, power_bounds,
    log_bounds, exact_floor_root, LOG_WITNESSES, ROOT_WITNESSES)
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
P=10**151
C,H=25237,Q(8793)
K,KL,L,CC=Q(159842),Q(662),Q(135),Q(26,25)
CM,CE=Q(26488),Q(26488)
DELTA=Q(12533,459000)
MID,HIGH,SWITCH=Q(10600),Q(15900),Q(3906)
EPS=Q(1,10**12)

def display(u,w,t):
    with localcontext() as ctx:
        ctx.prec=30
        def d(x):return Decimal(x.numerator)/Decimal(x.denominator)
        return str(((1258*d(u).ln()+377*d(w).ln()+252*d(w+EPS).ln()+252*d(t).ln())/1887).exp())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--moment',type=Path,required=True)
    ap.add_argument('--trace',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    if args.output.exists():raise SystemExit('Refusing to overwrite existing certificate')
    start=time.monotonic()
    result={'status':'RUNNING_FRACTIONAL_ENDPOINT_CANDIDATE','P':str(P),'C':C,'H_high':str(H),
      'fractional_reference':{'C_fractional':int(CM),'duplicated_reference_for_old_checker_only':int(CE),'R':3},
      'fixed_cover_switch':fraction(SWITCH),'cuts':[fraction(Q(128)),fraction(MID),fraction(HIGH)],
      'epsilon_upper':fraction(EPS),'scope':'Numerical endpoint comparisons; geometric reductions are proved in the accompanying article.',
      'gates':[],'endpoints':[],'bindings':[], 'root_acceptance_claimed':False,
      'script_sha256':file_sha(__file__),'helper_sha256':file_sha(HERE/'exact_arithmetic.py')}
    def persist():
        result.update(seconds=time.monotonic()-start,finished_utc=datetime.now(timezone.utc).isoformat(),
          log_witnesses=LOG_WITNESSES,root_witnesses=ROOT_WITNESSES)
        args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    def gate(name,condition,**detail):
        assert type(condition) is bool
        result['gates'].append({'name':name,'pass':condition,**detail})
        if not condition:raise AssertionError(name)
    def bind(rel):
        f=ROOT/rel
        result['bindings'].append({'path':rel,'sha256':file_sha(f),'bytes':f.stat().st_size})
        return f
    try:
        fractional=json.loads(args.moment.read_text())
        gate('fractional_moment_input',fractional['status']=='PASS_EXACT_FRACTIONAL_SLICE' and fractional['C_no3']==9525 and fractional['C_fixed_cover_r3_2']==26488 and fractional['all_q3_exponents_covered'])
        trace=json.loads(args.trace.read_text())
        gate('signed_trace_input',trace['status'].startswith('PASS_EDGEWORTH') and trace['frequency_cells']==15104)
        for row,target in zip(trace['cases'],[Q(17,50),Q(7,20)]):
            gate('trace_'+row['case'],Q(row['finite_maximum']['total_upper']['upper'])<target and Q(row['infinite_variance_tail']['upper'])<target)
        result['bindings']=[{'path':str(p),'sha256':file_sha(p)} for p in [args.moment,args.trace]]
        gate('fractional_budget_embedding',CM==CE==26488 and 1258+377+252==1887)
        sigma,eta,delta=Q(1,17),Q(1,18),Q(829,9180)
        r=(1+3*sigma+eta)/(1-3*sigma);kappa=3+3*r;alpha=(4+3*r)*delta
        gate('elimination_exponents',r==Q(377,252) and kappa==Q(629,84) and alpha==Q(591077,771120))
        gate('budget_weights',(2+2*r)/kappa==Q(1258,1887) and r/kappa==Q(377,1887) and 1/kappa==Q(252,1887))
        gate('conservative_final_powers',alpha<Q(174,227) and Q(125,252)>Q(41,84)
             and kappa*19068==142783 and Q(174,227)*19068==14616
             and (3+Q(41,84)+Q(174,227))*19068==81127)
        gate('epsilon_uniform_bound',Q(P,2)>10**150 and delta>Q(2,25)
             and Q(150)*Q(2,25)==12 and Q(1,3*10**12)<EPS)
        gate('delta_transport_exponent',0<DELTA<Q(7,256))
        gate('fixed_trace_transport',(Q(159,140)**256)>Q(100)**7
             and HIGH/Q(14000)==Q(159,140) and Q(10**153,P)==100 and MID==Q(2,3)*HIGH)
        left=int(H)**174960*2**71*P**469
        right=21718**174960*int(H+L)**540
        gate('ordinary_mass_capacity_cap',left>right,left=integer_digest(left),right=integer_digest(right))
        gate('capacity_increasing',324>1 and L>0 and Q(469,174960)>0)
        radicand=(P+1)**5;t=exact_floor_root(radicand,44);t+=(t**44<radicand)
        gate('ceil_root',(t-1)**44<radicand<=t**44,T=str(t))
        safe=36864*2541933**4*2541935**2*(t+1)**6
        cage=256*2541933**4*(t+1)**4
        gate('safe8_at_P',safe<P,left=str(safe),right=str(P))
        gate('cage_at_P',cage<P,left=str(cage),right=str(P))
        gate('uniform_safe_exponents',6*Q(5,44)-1==Q(-7,22) and 4*Q(5,44)-1==Q(-6,11))
        # Old safe8 -> old max-order -> upward-cover characteristic condition.
        gate('upward_characteristic_R3',648*3*2541935**2<36864*2541933**4)
        xh=K/KL;xl=K*CC/L
        gate('fixed_breakpoint_order',0<128<xh<xl<SWITCH<MID<HIGH<K)
        # We deliberately use weaker old U=W=H+low throughout the left branch;
        # no min-envelope or unrecorded moment hinge needs to be located.
        logp_lo,_=log_bounds(Q(P,2));_,log43540_hi=log_bounds(Q(43540))
        D128=log_bounds(Q(128))[0]+DELTA*logp_lo-log43540_hi
        gate('profile_log_convex',D128>3,D_lower=fraction(D128))
        points={'zero':Q(0),'128':Q(128),'horizontal':xh,'low_cap':xl,
                'switch':SWITCH,'mid':MID,'high':HIGH,'right':K}
        regions=[('old_trivial','old',['zero','128'],None),
                 ('old_profile','old',['128','horizontal','low_cap','switch'],None),
                 ('new_profile','new',['switch','mid'],None),
                 ('new_relaxed','new',['mid','high'],Q(7,20)),
                 ('new_strong','new',['high','right'],Q(17,50))]
        for region,cover,names,gamma in regions:
            for name in names:
                x=points[name];y=min(KL,K/x) if x else KL;low=min(L,CC*y)
                if cover=='old':u=w=H+low
                else:
                    powup=power_bounds(x,-2,15)[1]
                    u=CM*powup+low;w=CE*powup+low
                row={'region':region,'cover':cover,'name':name,'x':fraction(x),
                     'y':fraction(y),'low':fraction(low),'A_fractional_reference_first_copy':fraction(u),'A_fractional_reference_second_copy':fraction(w),'actual_standalone_mass_bound_claimed':False}
                if region.endswith('profile'):
                    dlo=log_bounds(x)[0]+DELTA*logp_lo-log43540_hi
                    gate('profile_D:'+region+':'+name,dlo>3)
                    gu=power_bounds(3/dlo,1,2)[1];tr=x*gu+y
                    row.update(D_lower=fraction(dlo),gamma_upper=fraction(gu))
                elif gamma is not None:tr=gamma*x+y;row['gamma']=fraction(gamma)
                else:tr=x+y
                q=u**1258*w**377*(w+EPS)**252*tr**252;rhs=C**1887
                row.update(trace_upper=fraction(tr),left=integer_digest(q.numerator),
                           right=integer_digest(rhs*q.denominator),F_display_only=display(u,w,tr),pass_exact=q<rhs)
                result['endpoints'].append(row);gate('geometric:'+region+':'+name,q<rhs)
        gate('all12_endpoints',len(result['endpoints'])==12)
        final_left=(3*C)**142783*(P+1)**14616
        final_right=2**81127*P**19068
        gate('final_exact_gate',final_left<final_right,left=integer_digest(final_left),right=integer_digest(final_right))
        gate('uniform_final_derivative',14616<19068)
        result['status']='PASS_FRACTIONAL_ENDPOINT_PENDING_INDEPENDENT_MATHEMATICAL_REVIEW'
        persist()
        print(json.dumps({k:result[k] for k in ['status','P','C','H_high','seconds']},indent=2))
        print(json.dumps({'gates':len(result['gates']),'endpoints':len(result['endpoints']),
                          'maximum_display_only':str(max(Decimal(z['F_display_only']) for z in result['endpoints']))}))
    except Exception as exc:
        result.update(status='FAILED_OR_ERROR_NO_ACCEPTANCE',error=repr(exc),traceback=traceback.format_exc());persist();raise

if __name__=='__main__':main()
