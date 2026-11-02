#!/usr/bin/env python3
"""Independent formulas and exact gates; imports no producer or producer helper.

Local moments use the integer-valued J first, then powers of log(q), and
variance uses E[J^2]-E[J]^2. Gaussian maxima use the original CDF expression
and fresh root brackets. Both arithmetic implementations share Arb/FLINT.
"""
from pathlib import Path
import os
from fractions import Fraction as F
import hashlib, json, math, sys, time
sys.set_int_max_str_digits(0)
ROOT = Path(os.environ['MARKOFF_WORKDIR'])
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
    p=(ROOT/dep['path'].split('/baragar_20261101_090837/',1)[1]); assert sha(p)==dep['sha256'] and p.stat().st_size==dep['bytes']
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
out={'status':'PASS_SELECTED_MOMENT_EULER_BUDGET_ONLY','finite_primes':len(prime_list),'finite_nonzero_exponents':count,'constant':43540,'script_sha256':sha(__file__)}
(ROOT/'root_logs/SELECTED_MOMENT_PUBLIC.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
