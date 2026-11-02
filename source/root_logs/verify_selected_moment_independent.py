#!/usr/bin/env python3
"""Independent formulas and exact gates; imports no producer or producer helper.

Local moments use the integer-valued J first, then powers of log(q), and
variance uses E[J^2]-E[J]^2. Gaussian maxima use the original CDF expression
and fresh root brackets. Both arithmetic implementations share Arb/FLINT.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math, sys, time
sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'agents/attachment_audit/d58/deps'))
from flint import arb, ctx
ctx.prec = 256
START = time.monotonic()
CERT = ROOT/'agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json'
c = json.loads(CERT.read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def aa(x):
    x=F(x); return arb(x.numerator)/x.denominator
def f(x): return F(int(x['numerator']),int(x['denominator']))
def box(x): return F(x['lower']),F(x['upper'])
def contained(v,saved):
    lo,hi=box(saved)
    assert aa(lo)<=v and v<=aa(hi)
def primes_below(n):
    mask=bytearray(b'\x01')*n
    mask[0:2]=b'\x00\x00'
    for p in range(2,math.isqrt(n-1)+1):
        if mask[p]: mask[p*p:n:p]=b'\x00'*len(range(p*p,n,p))
    return [p for p in range(2,n) if mask[p]]
def floor_root(n,k):
    low=0;high=1<<((n.bit_length()+k-1)//k)
    while high-low>1:
        m=(low+high)//2
        if m**k<=n:low=m
        else:high=m
    assert low**k<=n<(low+1)**k
    return low
def reciprocal_eighteenth(q):
    scale=1<<210
    z=floor_root((scale**18)//q,18)
    assert q*z**18<=scale**18<q*(z+1)**18
    return arb(aa(F(2*z+1,2*scale)),aa(F(1,2*scale)))
for dep in c['dependencies']:
    p=Path(dep['path']); assert sha(p)==dep['sha256'] and p.stat().st_size==dep['bytes']
e=c['euler']; assert e['prime_tail_from']==4099
prime_list=primes_below(4099)
assert [r['q'] for r in e['rows']]==prime_list and len(prime_list)==564
t=F(43,1000);u=F(61,10000);D=F(8747,1000);dp=F(63,1000)
delta=F(829,9180);P=3*10**153
prod=arb(1);count=0
for rec in e['rows']:
    q=rec['q'];lq=arb(q).log();r=reciprocal_eighteenth(q);cq=arb(q-1)/q
    inf=1+cq*r/(1-r)
    a0=rec['tail_first_exponent']
    assert [z['a'] for z in rec['values']]==list(range(1,a0))
    tail=(-aa(dp)*a0*lq).exp()*inf
    assert tail<1
    local=arb(1)
    for row in rec['values']:
        a=row['a'];weights=[arb(1)]+[cq*r**j for j in range(1,a+1)]
        B=sum(weights,arb(0))
        mu=sum((j*weights[j] for j in range(a+1)),arb(0))/B
        E2=sum((j*j*weights[j] for j in range(a+1)),arb(0))/B
        var=(E2-mu*mu)*lq*lq
        third=sum((weights[j]*abs(arb(j)-mu)**3 for j in range(a+1)),arb(0))/B*lq**3
        contained(B,row['B']);contained(mu*lq,row['mean'])
        contained(var,row['variance']);contained(third,row['third_absolute_centered'])
        aux=aa(t)*var-aa(u)*third
        penalty=aux.lower() if aux>0 else arb(0)
        value=(-aa(dp)*a*lq).exp()*B*(-penalty).exp()
        assert value<=aa(box(row['local_value_upper'])[1])
        upper=value.upper()
        if upper>local:local=upper
        count+=1
    assert local<=aa(box(rec['local_upper'])[1])
    prod*=local
assert count==1355==e['finite_nonzero_exponent_count']
assert prod<43540
assert (-aa(dp)*arb(4099).log()).exp()*(1+reciprocal_eighteenth(4099))<1
sizeD=arb(18000).log()+aa(delta-dp)*(arb(P)/2).log()-arb(43540).log()
assert sizeD>aa(D)
print('Independent local moments and all prime/exponent tails PASS',flush=True)
pi=arb.pi();root2=arb(2).sqrt()
def gaussian_max_at_variance(v):
    a=arb(v).sqrt()/18
    def cdf(w):return (1+(aa(w)/root2).erf())/2
    def mills(w):return (-aa(w)**2/2).exp()/((2*pi).sqrt()*cdf(w))
    lo=F(-8);hi=F(8)
    assert mills(lo)>a and mills(hi)<a
    for _ in range(60):
        mid=(lo+hi)/2
        if mills(mid)>a:lo=mid
        else:
            assert mills(mid)<a
            hi=mid
    assert mills(lo)>a and mills(hi)<a
    w=arb(aa((lo+hi)/2),aa((hi-lo)/2))
    z=a+w
    return (a*a/2-a*z).exp()*(1+(w/root2).erf())/2
vd=c['variance_domain']; assert vd['variance_boxes']==997 and D/t>203
assert [x['variance_interval'] for x in vd['rows']]==[[v,v+1] for v in range(203,1200)]
max_bound=arb(0)
for rec in vd['rows']:
    v0,v1=rec['variance_interval'];H=gaussian_max_at_variance(v0)
    critical=3*D/t;arg=min(F(v1),max(F(v0),critical))
    assert F(rec['error_argmax'])==arg
    error=aa(F(28,25)/u)*(aa(t*arg-D))/(aa(arg)**3).sqrt()
    if error<0:error=arb(0)
    upper=H.upper()+error.upper()
    assert upper<aa(F(97,200))
    if upper>max_bound:max_bound=upper
tail=(18/(2*pi).sqrt()+aa(F(28,25)*t/u))/arb(1200).sqrt()
assert tail<aa(F(97,200))
new_sum=F(97,200)+F(103,200)*F(159842,18000**2+159842)
assert new_sum==f(vd['sum_ratio_upper']) and new_sum<F(243,500)
print('Fresh normal roots,997 continuous boxes and infinite tail PASS',flush=True)
def exp_upper(z):return sum(((-z)**j/math.factorial(j) for j in range(17)),F(0))
def log_lower(z):
    assert z>=1
    a=(z-1)/(z+1)
    return 2*sum((a**(2*j+1)/F(2*j+1) for j in range(45)),F(0))
def atan_bracket(z):
    val=sum(((-1)**j*z**(2*j+1)/F(2*j+1) for j in range(41)),F(0))
    return val-z**83/83,val
g=c['old_probability_gate'];old=ROOT/'agents/analytic_window_round4/SIZE_TRACE_CERTIFICATE.json'
accepted=json.loads((ROOT/'root_logs/SIZE_TRACE_ROOT_REVIEW.json').read_text())
assert accepted['status']=='PASS_INDEPENDENT_EXACT_SIZE_TRACE_AND_ENDPOINT'
assert sha(old)==accepted['certificate_sha256']==g['prior_size_certificate_sha256']
assert 40300**18360*2**281<183**18360*P**281
al,au=atan_bracket(F(1,5));bl,bu=atan_bracket(F(1,239));piL=16*al-4*bu
ll=5*log_lower(F(2))+log_lower(F(10000,32*183));q=F(1411,10000)
qs=F(1,4)/(piL*ll);assert qs==f(g['q_squared_upper']) and qs<q*q
su=exp_upper(F(5,24));assert 7*q<=1<8*q
gamma=q*sum((su**j for j in range(7)),F(0))+(1-7*q)*su**7
assert gamma==f(g['gamma_upper'])
ratio=gamma+(1-gamma)*F(159842,10000**2+159842)
assert ratio==f(g['sum_ratio_upper']) and ratio<F(29,50)
end=c['endpoint'];assert int(end['P'])==P and end['C']==30130 and end['ordinary_high_cap']==8840
T=int(end['prerequisites']['T']);assert (T-1)**44<(P+1)**5<=T**44
assert 36864*2541933**4*2541935**2*T**6<P
assert 256*2541933**4*T**4<P
assert 21718**201960*2541935**660<8840**201960*2**119*P**541
assert (729*(2*181265+3))**900*27**420*14556**6660*(P+1)**673<(16*P)**900
H=F(8840);L=F(301);Cm=F(34195);CE=F(28924);cc=F(11,10);K=F(159842);KL=F(1826);C=30130
expected=[('low',x) for x in ['left_limit','horizontal_to_hyperbola','low_count_crossover','high_mass_crossover','split_x_10000']]
expected += [('middle','split_x_10000'),('middle','split_x_18000'),('high','split_x_18000'),('high','high_count_crossover'),('high','right_edge')]
pts=end['point_bounds'];assert [(x['region'],x['point']) for x in pts]==expected
seen={}
for rec in pts:
    x,y,U,Vm=[f(rec[k]) for k in ['x_upper','y_upper','U_upper','V_mass_upper']]
    label=rec['point'];fac={'low':F(1),'middle':F(29,50),'high':F(243,500)}[rec['region']]
    assert F(rec['trace_ratio_used'])==fac
    if label=='left_limit':assert (x,y,U,Vm)==(0,KL,H+L,H+L)
    elif label=='horizontal_to_hyperbola':assert (x,y,U,Vm)==(K/KL,KL,H+L,H+L)
    elif label=='low_count_crossover':assert (x,y,U,Vm)==(K/(L/cc),L/cc,H+L,H+L)
    elif label in ['high_mass_crossover','high_count_crossover']:
        lower=K/y;target=(CE/H if label=='high_mass_crossover' else Cm/H)**15
        assert lower**2<=target<=x**2 and x-lower==F(1,1<<80)
        assert U==H+cc*y
        assert Vm==(H+cc*y if label=='high_mass_crossover' else H*CE/Cm+cc*y)
        if label=='high_mass_crossover':assert x<10000
        else:assert lower>18000
    elif label.startswith('split_x_'):
        A=F(label.split('_')[-1]);assert A in [10000,18000] and x==A and y==K/A
        assert U==H+cc*y
        power=(Vm-cc*y)/CE
        assert power>0 and power**15*A**2>=1 and CE*power<H
    elif label=='right_edge':
        assert x==K and y==1
        power=(U-cc)/Cm;assert power>0 and power**15*K**2>=1 and Vm==CE*power+cc
    else:raise AssertionError(label)
    value=U**1510*Vm**377*((x+y)*fac)**252
    hi=f(rec['F_upper']);assert value<=hi**1887 and hi<C
    if label in seen:assert seen[label]==(x,y,U,Vm)
    seen[label]=(x,y,U,Vm)
left=(3*C)**142783*(P+1)**14616;right=2**81127*P**19068
assert left<right and [left.bit_length(),right.bit_length()]==end['endpoint_comparison_bits']
out={'status':'PASS_INDEPENDENT_SELECTED_MOMENT_FULL_ARITHMETIC_AND_ENDPOINT',
     'certificate_sha256':sha(CERT),'verifier_sha256':sha(__file__),'precision_bits':ctx.prec,
     'fresh_prime_factors':len(prime_list),'fresh_nonzero_exponents':count,
     'fresh_normal_root_brackets':997,'continuous_variance_boxes':997,'infinite_variance_tail_checked':True,
     'new_euler_product_display_only':str(prod),'maximum_variance_bound_display_only':str(max_bound),
     'endpoint_points':len(pts),'P':str(P),'endpoint_binary_digits':[left.bit_length(),right.bit_length()],
     'seconds':time.monotonic()-START,
     'scope':'No producer/helper imports. Different moment and Gaussian formulas at256bits; all finite terms and all tails rebuilt. Exact Fraction endpoint checks. Shared Arb/FLINT is disclosed; external nonidentical Berry-Esseen C0.56 and inherited Markoff geometry remain separately reviewed theorem dependencies.'}
(ROOT/'root_logs/SELECTED_MOMENT_ROOT_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
