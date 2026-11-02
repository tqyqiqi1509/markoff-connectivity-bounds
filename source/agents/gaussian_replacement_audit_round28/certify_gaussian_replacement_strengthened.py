#!/usr/bin/env python3
"""One new continuous-variance run after the centered 4b symmetrization proof.

Reuses this agent's immutable exact frequency table and Gaussian helper.
This reuse is explicit; independent replay/acceptance is the root's task.
No numerical producer from the inherited selected-moment work is imported.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util, json, math, time

HERE=Path(__file__).resolve().parent
BASE=HERE/'certify_gaussian_replacement.py'
spec=importlib.util.spec_from_file_location('gaussian_replacement_base',BASE)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
arb,aa,pack=m.arb,m.aa,m.pack
GAMMA=F(21,50)

def one_sided_normal(v):
    a=aa(v).sqrt()/18
    def mills(w):
        z=aa(w)
        return (-z*z/2).exp()/((2*m.PI).sqrt()*((-z/m.SQRT2).erfc()/2))
    lo=F(-10);hi=F(10)
    assert mills(lo)>a and mills(hi)<a
    for _ in range(45):
        mid=(lo+hi)/2
        mm=mills(mid)
        if mm>a: lo=mid
        elif mm<a: hi=mid
        else: raise ArithmeticError('Unseparated BE normal maximum root')
    assert mills(lo)>a and mills(hi)<a
    w=m.box(lo,hi)
    value=(-(a+w)**2/2).exp()/((2*m.PI).sqrt()*a)
    return value.upper(),lo,hi

def main():
    started=time.monotonic()
    tablepath=HERE/'EULER_FREQUENCY_CERTIFICATE.json'
    table=json.loads(tablepath.read_text())
    assert table['parameters']==m.parameters()
    assert table['status']=='PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES'
    assert table['method']=='Directed midpoint finite sums plus rigorous per-candidate global derivative times half width.'
    # The prior target was run first and remains a separately bound result.
    previouspath=HERE/'GAUSSIAN_REPLACEMENT_CERTIFICATE.json'
    previous=json.loads(previouspath.read_text())
    assert previous['status']=='PASS_NEW_GAUSSIAN_REPLACEMENT_UNIFORM_TRACE_PENDING_ROOT_REVIEW'
    assert previous['euler_table_sha256']==m.digest(tablepath)
    assert previous['producer_sha256']==m.digest(BASE)
    gatepath=HERE/'NEW_P_BUDGET_GATE.json'
    gate=json.loads(gatepath.read_text())
    assert gate['status']=='PASS_INHERITED_K_AND_NEW_P_D_GATE'
    assert F(gate['D_lower'])==m.D
    assert F(gate['D_expression_lower'])>m.D
    grids=m.frequency_grid(); values=[aa(F(x)) for x in table['frequency_upper_values']]
    assert len(grids)==len(values)==9984
    cells=[]
    for (lo,hi),aval in zip(grids,values):
        weight=(aa(m.ETA+m.BETA)/(m.PI*((aa(m.ETA)**2+aa(lo)**2)*(aa(m.BETA)**2+aa(lo)**2)).sqrt())).upper()
        cells.append((lo,hi,aval,weight*aa(hi-lo)))
    lowcells=cells[:2048];highcells=cells[2048:]
    highA=sum((aval*mass for _,_,aval,mass in highcells),arb(0)).upper()
    floor=m.D/m.LAM;start=math.floor(floor);stop=1600
    frequency_tail=((aa(m.ETA+m.BETA)/(8*m.PI))*(1+(-aa(floor)*32).exp())).upper()
    critical=3*m.D/m.LAM
    maxupper=arb(0);worst=None;rows=[];used_BE=0
    for vl in range(start,stop):
        v0=max(F(vl),floor);v1=F(vl+1)
        if v0>=v1: continue
        Tmax=aa((m.LAM*v1-m.D)/m.NU)
        C,lroot,rroot=m.normal_peak(v0)
        total=arb(0)
        for lo,hi,aval,mass in lowcells:
            l=aa(lo);u=aa(hi)
            gaussian=(-aa(v0)*l*l/2).exp().upper()
            cap=(aval+gaussian).upper()
            pref=Tmax*(u**3/6+aa(m.C)*u**4/8)
            exponent=-aa(v0)*l*l/2+Tmax*u**3/3+aa(m.C)**2*u*u/2
            err=(pref*exponent.exp()).upper().min(cap).upper()
            total+=mass*err
        highG=sum((mass*(-aa(v0)*aa(lo)**2/2).exp().upper() for lo,_,_,mass in highcells),arb(0)).upper()
        finite_error=total.upper()+highA+highG
        replacement=(C+finite_error+frequency_tail).upper()
        H,wlo,whi=one_sided_normal(v0)
        arg=min(v1,max(v0,critical))
        beerror=(aa(F(28,25)/m.NU)*aa(m.LAM*arg-m.D)/(aa(arg)**3).sqrt()).upper()
        assert beerror>=0
        be=(H+beerror).upper()
        full=replacement.min(be).upper()
        chosen='BE' if be<replacement else 'Gaussian_replacement'
        if chosen=='BE': used_BE+=1
        if not full<aa(GAMMA):
            failure={'status':'FAILED_STRENGTHENED_TARGET','V_interval':[str(v0),str(v1)],
                     'replacement':pack(replacement),'BE':pack(be),'total':pack(full)}
            m.save('STRENGTHENED_VARIANCE_TARGET_FAILURE.json',failure)
            raise ArithmeticError(json.dumps(failure))
        if full>maxupper: maxupper=full;worst=[str(v0),str(v1)]
        rows.append({'V_interval':[str(v0),str(v1)],'T_upper':pack(Tmax),
            'Gaussian_root_interval':[str(lroot),str(rroot)],'Gaussian_peak_upper':pack(C),
            'finite_error_upper':pack(finite_error),'replacement_upper':pack(replacement),
            'BE_normal_root_interval':[str(wlo),str(whi)],'BE_normal_H_upper':pack(H),
            'BE_error_argmax':str(arg),'BE_error_upper':pack(beerror),'BE_upper':pack(be),
            'chosen_complete_bound':chosen,'total_upper':pack(full)})
        if vl%100==0:
            print('Strengthened variance boxes through',vl,'max',str(maxupper),'elapsed',round(time.monotonic()-started,2),flush=True)
    infinite_tail=m.be_tail(stop)
    assert infinite_tail<aa(GAMMA)
    global_upper=maxupper.max(infinite_tail.upper()).upper()
    exact_upper=F(str(global_upper.fmpq()))
    rounded_gamma=F((1000*exact_upper.numerator)//exact_upper.denominator+1,1000)
    assert global_upper<aa(rounded_gamma) and rounded_gamma<=GAMMA
    params=m.parameters();params['gamma_target']=str(GAMMA);params['cubic_exponent_coefficient']='1/3'
    dependencies=[BASE,Path(__file__),tablepath,previouspath,gatepath,
        HERE/'INDEPENDENT_GAUSSIAN_REPLACEMENT_PROOF.md',HERE/'CENTERED_SYMMETRIZATION_STRENGTHENING.md',
        m.ROOT/'root_logs/SELECTED_MOMENT_ROOT_REVIEW.json',
        m.ROOT/'agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json']
    out={'status':'PASS_STRENGTHENED_GAUSSIAN_REPLACEMENT_ALL_VARIANCES_PENDING_ROOT_REVIEW',
         'parameters':params,'precision_bits':m.ctx.prec,
         'variance_floor':str(floor),'variance_boxes':len(rows),'variance_box_width':1,
         'maximum_box':worst,'maximum_upper':pack(maxupper),'boxes_using_BE':used_BE,
         'gamma_certified_after_exact_rounding':str(rounded_gamma),
         'rounding_rule':'Smallest strict upper grid value with denominator 1000, applied once to the completed certificate maximum including the infinite variance tail.',
         'finite_high_frequency_A_upper':pack(highA),'frequency_infinite_tail_upper':pack(frequency_tail),
         'infinite_variance_tail_from':stop,'infinite_variance_tail_upper':pack(infinite_tail),
         'variance_rows':rows,'all_variances_covered':True,'all_frequencies_covered':True,
         'all_primes_exponents_covered_by_table':True,'new_external_theorems':0,
         'inherited_external_theorem':'Nonidentical independent Berry-Esseen C=0.56, already source-audited and root accepted.',
         'dependencies':[{'path':str(p),'sha256':m.digest(p),'bytes':p.stat().st_size} for p in dependencies],
         'seconds':time.monotonic()-started,
         'scope':'Fixed P=10^153 sufficient for the probability estimate; no geometric endpoint is asserted here.'}
    m.save('STRENGTHENED_GAUSSIAN_REPLACEMENT_CERTIFICATE.json',out)
    print(json.dumps({k:out[k] for k in ['status','maximum_box','maximum_upper','gamma_certified_after_exact_rounding','boxes_using_BE','seconds']},indent=2),flush=True)

if __name__=='__main__':main()
