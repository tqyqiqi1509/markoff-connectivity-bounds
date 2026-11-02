#!/usr/bin/env python3
"""Size-preserving enhancement; exact arithmetic, no new spectral nodes.

Imports shared elementary routines from the preceding producer. This is a
new certificate, NOT an independent verification of that producer.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,sys,time,hashlib

sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
PREVIOUS=HERE/'certify_trace_concentration.py'
spec=importlib.util.spec_from_file_location('trace_preceding_producer',PREVIOUS)
OLD=importlib.util.module_from_spec(spec);spec.loader.exec_module(OLD)
HLP=OLD.HLP;S=OLD.S;pack=OLD.pack;ceildiv=OLD.ceildiv;upmul=OLD.upmul

def size_euler():
    elo,ehi=OLD.exp_pos_bounds(F(15,4))
    assert 42<elo<ehi<43
    rows=[];product=S
    for q in HLP.primes_before(1019):
        wl,wu=HLP.prime_power_box(q,3,40)
        rl,ru=HLP.prime_power_box(q,1,18)
        k=1
        while q**k<43:k+=1
        assert q**(k-1)<=42 and q**k>=43
        binf=S+ceildiv((q-1)*ru*S,q*(S-ru))
        wp=rp_l=rp_u=S;sum_l=sum_u=0;local=S
        values=[{'a':0,'lambda_lower_scaled':0,'value_upper_scaled':S}]
        a=0;num_u=S
        while True:
            a+=1;wp=upmul(wp,wu)
            rp_l=rp_l*rl//S;rp_u=upmul(rp_u,ru)
            sum_l+=rp_l;sum_u+=rp_u
            if wp*binf<=S*S:break
            bl=S+(q-1)*sum_l//q
            bu=S+ceildiv((q-1)*sum_u,q)
            if a<k:num_u=bu;lam=0
            else:
                qh_u=ceildiv(num_u*S,bl)
                lam=max(0,min(S//2,S-qh_u))
            ex=OLD.exp_minus_upper(F(lam,4*S))
            exu=ceildiv(ex.numerator*S,ex.denominator)
            val=upmul(upmul(wp,bu),exu);local=max(local,val)
            values.append({'a':a,'B_lower_scaled':bl,'B_upper_scaled':bu,
                           'first_k_numerator_upper_scaled':num_u,
                           'lambda_lower_scaled':lam,'exp_upper_scaled':exu,
                           'value_upper_scaled':val})
            assert a<1000
        product=upmul(product,local)
        rows.append({'q':q,'k_atoms':k,'w_box':[wl,wu],'r_box':[rl,ru],
                     'local_upper_scaled':local,'finite_values':values,
                     'tail_first_exponent':a,'tail_w_power_upper_scaled':wp,
                     'B_infinite_upper_scaled':binf})
    wl,wu=HLP.prime_power_box(1019,3,40)
    rl,ru=HLP.prime_power_box(1019,1,18)
    assert wu*(S+ru)<S*S
    assert product<40300*S
    delta=F(829,9180)-F(3,40)
    assert delta==F(281,18360)
    assert 40300**18360*2**281<180**18360*(10**154)**281
    return {'status':'EXACT_SIZE_EULER_PASSED','scale':S,'h':'15/4','t':'1/4',
            'eta':'1/18','delta_prime':'3/40','delta_reference':'829/9180',
            'delta_difference':str(delta),'integer_cap':40300,
            'K_upper':pack(F(product,S)),'K_display_only':float(F(product,S)),
            'exp_h_bracket':{'lower':pack(elo),'upper':pack(ehi),'integer_bracket':[42,43]},
            'rows':rows,'infinite_prime_tail':{'first_prime':1019,'w_box':[wl,wu],
            'r_box':[rl,ru],'criterion':'w_upper*(1+r_upper)<1; decreases with q'},
            'size_gate':{'P0':str(10**154),'effective_cap':180,
            'left':'40300^18360*2^281','right':'180^18360*(10^154)^281','strict_comparison':True}}

def size_factor():
    al,au=OLD.atan_bounds(F(1,5));bl,bu=OLD.atan_bounds(F(1,239))
    pi_l=16*al-4*bu;pi_u=16*au-4*bl
    assert 3<pi_l<pi_u<F(22,7)
    expu=OLD.exp_minus_upper(F(5,24))
    log_l=6*OLD.log_lower(F(2))+OLD.log_lower(F(17000,64*180))
    bound=F(1,4)/(pi_l*(1-expu)**2*log_l)
    rho=F(88,125);ratio=F(141,200)
    assert bound<rho*rho
    assert rho+(1-rho)*F(159842,17000**2+159842)<ratio
    return {'status':'EXACT_SIZE_TRACE_FACTOR_PASSED','pi_lower':pack(pi_l),'pi_upper':pack(pi_u),
            'exp_minus_5_over_24_upper':pack(expu),'log_17000_over_180_lower':pack(log_l),
            'rho_squared_bound':pack(bound),'rho':str(rho),'A':17000,
            'sum_ratio':str(ratio),'rho_bound_display_only':float(bound)**.5}

def size_endpoint():
    P=15*10**153;H=F(8802);L=F(301);Cm=F(34195);CE=F(28924)
    K=F(159842);KL=F(1826);c=F(11,10);A=F(17000);ratio=F(141,200);C=31700
    ir=HLP.integer_root((P+1)**5,44);T=ir+int(ir**44!=(P+1)**5)
    assert (T-1)**44<(P+1)**5<=T**44
    safe=36864*2541933**4*2541935**2*T**6;cage=256*2541933**4*T**4
    assert safe<P and cage<P
    assert 21718**201960*2541935**660<int(H)**201960*2**119*P**541
    p63=10**154
    assert P>=p63
    assert (729*(2*181265+3))**900*27**420*14556**6660*(p63+1)**673<(16*p63)**900
    eL,eU=HLP.power_box(CE/H,15,2);mL,mU=HLP.power_box(Cm/H,15,2)
    xh=K/KL;yl=L/c;xl=K/yl
    assert 0<xh<xl<eL<=eU<A<mL<=mU<K
    _,rp=HLP.power_box(K,-2,15)
    pts=[('left_limit',F(0),KL,H+L,H+L),
         ('horizontal_to_hyperbola',xh,KL,H+L,H+L),
         ('low_count_crossover',xl,yl,H+L,H+L),
         ('high_mass_crossover',eU,K/eL,H+c*K/eL,H+c*K/eL),
         ('high_count_crossover',mU,K/mL,H+c*K/mL,H*CE/Cm+c*K/mL),
         ('right_edge',K,F(1),Cm*rp+c,CE*rp+c)]
    _,ap=HLP.power_box(A,-2,15)
    at=('split_x_17000',A,K/A,min(H,Cm*ap)+c*K/A,min(H,CE*ap)+c*K/A)
    rows=[]
    for name,items,factor in [('low',pts[:4]+[at],F(1)),('high',[at]+pts[4:],ratio)]:
        for label,x,y,U,V in items:
            base=U**1510*V**377*((x+y)*factor)**252
            assert base<C**1887
            lo,hi=HLP.power_box(base,1,1887)
            rows.append({'region':name,'point':label,'x_upper':pack(x),'y_upper':pack(y),
                         'U_upper':pack(U),'V_upper':pack(V),'trace_ratio_used':str(factor),
                         'F_upper':pack(hi),'F_display_only':float(hi)})
    left=(3*C)**142783*(P+1)**14616;right=2**81127*P**19068
    assert left<right
    return {'status':'EXACT_SIZE_ENDPOINT_PASSED_CONDITIONAL_ON_GENERAL_PROOF','P':str(P),
            'C':C,'ordinary_high_cap':int(H),'point_bounds':rows,
            'gates':{'safe8':True,'cage':True,'ordinary_transfer':True,
                     'other_63_at_1e154':True,'two_region_full_domain':True,
                     'endpoint_integer_comparison':True},
            'endpoint':{'left_formula':f'(3*{C})^142783*(15*10^153+1)^14616',
                        'right_formula':'2^81127*(15*10^153)^19068',
                        'left_bits':left.bit_length(),'right_bits':right.bit_length()},
            'prerequisites':{'T':str(T),'safe8_left':str(safe),'cage_left':str(cage)}}

def main():
    started=time.monotonic()
    result={'status':'PENDING'}
    result['euler']=size_euler();result['trace_factor']=size_factor();result['endpoint']=size_endpoint()
    result['status']='EXACT_SIZE_CERTIFICATE_COMPLETE_PENDING_ROOT_REVIEW'
    result['scope']='Uniform sufficient bound, conditional on accepted preceding geometry and spectral inputs plus KR/Kanter and the new size proof. Not the full conjecture.'
    deps=[Path(__file__),PREVIOUS,OLD.HELPER,HERE/'FINITE_SIZE_TRACE_ENHANCEMENT_PROOF.md',
          HERE/'TRUE_TRUNCATED_TRACE_CONCENTRATION_PROOF.md',
          HERE.parent/'analytic_mass_round3/MASS_MOMENT_CERTIFICATE.json',
          HERE.parent/'joint_analytic/HALF_MODULUS_ENDPOINT_CERTIFICATE.json']
    result['dependencies']=[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in deps]
    result['seconds']=time.monotonic()-started
    (HERE/'SIZE_TRACE_CERTIFICATE.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'K_upper_display_only':result['euler']['K_display_only'],
                     'effective_cap':180,'rho_bound_display_only':result['trace_factor']['rho_bound_display_only'],
                     'endpoint':{k:v for k,v in result['endpoint'].items() if k not in ['point_bounds','prerequisites']},
                     'points':[{k:r[k] for k in ['region','point','F_display_only']} for r in result['endpoint']['point_bounds']],
                     'seconds':result['seconds']},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
