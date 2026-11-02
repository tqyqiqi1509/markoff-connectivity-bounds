#!/usr/bin/env python3
"""Readback audit, separate rational recurrence, and tilted trace Euler bound.

This program does not import the integration program.  Its recurrence uses
Fractions with explicit rational upward rounding, not the integrator's Arb
recurrence.  Arb at 384 bits only certifies roots and local analytic data.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path
import struct
import sys
import time

import flint
from flint import arb, ctx

sys.set_int_max_str_digits(0)
ctx.prec = 384


def fraction(d):
    m, e = int(d["mantissa"]), d["exponent"]
    return F(m*(1<<e)) if e >= 0 else F(m,1<<(-e))


def dyad(x):
    assert x.is_exact()
    m,e = x.man_exp()
    return {"mantissa":str(m),"exponent":int(e)}


def ball(x):
    return arb(x.numerator)/x.denominator


def ceil_fraction(x):
    return -(-x.numerator//x.denominator)


def up(x,bits=96):
    return F(ceil_fraction(x*(1<<bits)),1<<bits)


def small_primes(limit):
    return [n for n in range(2,limit+1)
            if all(n%d for d in range(2,isqrt(n)+1))]


def ipow(q,t):
    return arb(q).root(t.denominator)**(-t.numerator)


def max_count(ds, num, den):
    small = [v for v in ds if den*v <= num]
    return sum(not any(u>v and u%v==0 for u in small) for v in small)


def structural_checks(limit=4096):
    ds = [[] for _ in range(limit+1)]
    for v in range(1,limit+1):
        for n in range(v,limit+1,v):
            ds[n].append(v)
    tests=0
    for n in range(2,limit+1):
        q=ds[n][1]
        qa, np=q,n//q
        while np%q==0:
            qa*=q
            np//=q
        cuts=[1,2*n+1]+[2*v for v in ds[n]]
        cuts += [a+b for a,b in zip(ds[n],ds[n][1:])]
        for num in cuts:
            lhs=max_count(ds[n],num,2)
            rhs=sum(num<2*qa*v and 2*v<=num for v in ds[np])
            rhs+=max_count(ds[np],num,2*qa)
            assert lhs==rhs,(n,num,2,lhs,rhs)
            tests+=1
    # Deliberate removal of the minimum-prime hypothesis fails.
    lhs=max_count(ds[6],2,1)
    rhs=sum(2<3*v and v<=2 for v in ds[2])+max_count(ds[2],2,3)
    assert lhs==1 and rhs==2
    return {"integer_n_range":[2,limit],"exact_cutoff_checks":tests,
            "all_positive_real_cutoff_cells_per_n":True,
            "counterexample_without_minimum_prime":{"n":6,"q":3,"a":1,"d":2,"lhs":lhs,"rhs":rhs}}


def validate_local(data):
    eps,sig=F(data["epsilon"]),F(data["sigma"])
    rows=data["parameters"]["active_primes"]
    expected=small_primes(data["parameters"]["first_omitted_prime"]-1)
    assert [r["q"] for r in rows]==expected
    for row in rows:
        q,E=row["q"],row["first_discarded_exponent"]
        w,r=ipow(q,eps),ipow(q,sig)
        for key,value in [("w",w),("r",r)]:
            assert value.lower() >= ball(fraction(row[key+"_lower"]))
            assert value.upper() <= ball(fraction(row[key+"_upper"]))
        assert (w**E/(1-r)).upper() <= ball(fraction(row["discarded_tail_upper"])) < 1
        rp,S,D=arb(1),arb(1),arb(0)
        for e in range(1,E):
            rp*=r
            S+=rp
            D+=e*rp
            assert (w**e*S).upper() <= ball(fraction(row["K_upper"]))
            assert (w**e*D).upper() <= ball(fraction(row["L_upper"]))
    q0=data["parameters"]["first_omitted_prime"]
    assert (ipow(q0,eps)*(1+ipow(q0,sig))).upper() < 1
    return len(rows)


def readback(data,path):
    info=data["integration"]
    n=len(info["tails"])
    codec=struct.Struct("<"+"Q"*n)
    assert codec.size==info["row_bytes"]
    scale=1<<info["node_scale_bits"]
    sums=[0]*n
    block_sums=[0]*n
    digest=hashlib.sha256()
    block_index=0
    raw=path.parent/info["node_file"]
    assert hashlib.sha256(raw.read_bytes()).hexdigest()==info["node_gzip_sha256"]
    with gzip.open(raw,"rb") as f:
        for i in range(info["panels"]):
            b=f.read(codec.size)
            assert len(b)==codec.size
            digest.update(b)
            values=codec.unpack(b)
            assert values[-1]==scale and all(v>=scale for v in values)
            assert all(x>=y for x,y in zip(values,values[1:]))
            for j,v in enumerate(values):
                sums[j]+=v
                block_sums[j]+=v
            if i==info["blocks"][block_index]["last_node"]:
                assert block_sums==info["blocks"][block_index]["tail_sums"]
                block_sums=[0]*n
                block_index+=1
        assert f.read(1)==b""
    assert block_index==len(info["blocks"])
    assert digest.hexdigest()==info["node_raw_sha256"]
    assert sums==[int(t["node_sum_integer"]) for t in info["tails"]]
    correction,K=arb(1),arb(1)
    params=data["parameters"]["active_primes"]
    for j in range(n-2,-1,-1):
        L=ball(fraction(params[j]["L_upper"]))
        correction*=1+L*arb.pi()/(2*info["panels"])
        K*=ball(fraction(params[j]["K_upper"]))
        tail=info["tails"][j]
        assert correction.upper()<=ball(fraction(tail["correction_upper"]))
        assert K.upper()<=ball(fraction(tail["K_upper"]))
        integral=arb(sums[j])/(info["panels"]*scale)*correction
        legitimate=min(integral.upper(),K.upper())
        assert legitimate<=ball(fraction(tail["W_upper"]))
    return {"all_nodes":info["panels"],"all_tails":n,
            "all_saved_integer_values":info["panels"]*n,
            "blocks":len(info["blocks"]),"sha256":digest.hexdigest(),
            "all_values_and_sums_verified":True}


def rational_recursion(data):
    rows=data["parameters"]["active_primes"]
    tails=data["integration"]["tails"]
    n=len(rows)
    K=[fraction(r["K_upper"]) for r in rows]
    W=[fraction(t["W_upper"]) for t in tails]
    Ktails=[F(1)]*(n+1)
    for j in range(n-1,-1,-1):
        Ktails[j]=up(K[j]*Ktails[j+1])
    C=F(1)
    exact_two=None
    results=[]
    for j in range(n-1,-1,-1):
        q=rows[j]["q"]
        w,r=(fraction(rows[j][x+"_upper"]) for x in ("w","r"))
        a=1
        while w**a*(Ktails[j+1]+C)>C:
            a+=1
        A=a
        best,best_a=C,0
        for a in range(1,A):
            J=Ktails[j+1]
            E=F(1)
            for k in range(j+1,n):
                Q=rows[k]["q"]
                h=1
                while Q**h < q**a:
                    h+=1
                J=min(J,up(E*h*W[k]))
                E=up(E*K[k])
            value=up(w**a*(J+r**a*C))
            if j==0 and a==1:
                exact_two=value
            if value>best:
                best,best_a=value,a
        C=best
        results.append({"q":q,"C_numerator":str(C.numerator),"C_denominator":str(C.denominator),
                        "first_omitted_a":A,"best_a":best_a})
    joint=C+exact_two
    integer=ceil_fraction(joint)
    assert integer==data["recursion"]["C_joint_integer_upper"]
    return integer,{"arithmetic":"exact rational, upward rounding to multiples of 2^-96",
        "C_joint_integer_upper":integer,
        "C_joint_numerator":str(joint.numerator),"C_joint_denominator":str(joint.denominator),
        "C_joint_display_only":ball(joint).str(40),"rows":results}


def tilted_trace_constant(delta=F(9,100),eta=F(1,18)):
    k=(delta.denominator+delta.numerator-1)//delta.numerator
    limit=1<<k
    assert limit**delta.numerator>=2**delta.denominator
    product=arb(1)
    records=[]
    active=[]
    for q in small_primes(limit):
        w,r=ipow(q,delta),ipow(q,eta)
        c=arb(q-1)/q
        e,B,power=0,arb(1),r
        while True:
            nextB=B+c*power
            ratio=w*nextB/B
            if ratio>1:
                B=nextB
                e+=1
                power*=r
            elif ratio<1:
                break
            else:
                raise ArithmeticError("undecided local maximum")
        factor=w**e*B
        if e:
            active.append(q)
            previousB=B-c*r**e
            assert (w*B/previousB).lower()>1
        assert ratio.upper()<1
        product*=factor
        records.append({"q":q,"maximizer_e":e,"ratio_after_upper":dyad(ratio.upper()),
                        "factor_lower":dyad(factor.lower()),"factor_upper":dyad(factor.upper())})
    integer=int(product.upper().ceil().fmpz())
    assert integer<=181265
    return integer,{"delta":str(delta),"eta":str(eta),"prime_limit":limit,
        "active_primes":active,"local_factors":records,
        "K_B_lower":dyad(product.lower()),"K_B_upper":dyad(product.upper()),
        "K_B_display_only":product.str(45),"K_B_integer_upper":integer}


def augmented_endpoint(C,K,power):
    P=10**power
    lhs=(216*(2*K+3))**900 * 27**420 * C**6660 * (P+1)**673
    rhs=(2*P)**900
    T=int((arb(P+1).root(44)**5).upper().ceil().fmpz())
    assert (T-1)**44 < (P+1)**5 <= T**44
    safe=36864*2541933**4*2541935**2*T**6
    cage=256*2541933**4*T**4
    ratio=arb(216)*(arb(K)+arb(3)/2)*arb(27)**(arb(7)/15)*arb(C)**(arb(37)/5)*arb(P+1)**(arb(673)/900)/P
    return {"power":power,"C_joint_integer_upper":C,"K_B_integer_upper":K,
        "checks":{"main":lhs<rhs,"safe8":safe<P,"cage":cage<P,"Arb_main":ratio<1},
        "main_ratio_upper":dyad(ratio.upper()),"main_ratio_display_only":ratio.str(45),
        "T":str(T),"safe8_left":str(safe),"cage_left":str(cage),
        "main_left_sha256":hashlib.sha256(str(lhs).encode()).hexdigest(),
        "main_right_sha256":hashlib.sha256(str(rhs).encode()).hexdigest()}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--structural-limit",type=int,default=4096)
    args=ap.parse_args()
    path=Path(args.input)
    data=json.loads(path.read_text())
    start=time.monotonic()
    source=path.parent/"independent_interval_arb.py"
    assert hashlib.sha256(source.read_bytes()).hexdigest()==data["script_sha256"]
    count=validate_local(data)
    saved=readback(data,path)
    integer,recur=rational_recursion(data)
    structure=structural_checks(args.structural_limit)
    K,KB=tilted_trace_constant()
    ep=augmented_endpoint(integer,K,176)
    ep_lower=augmented_endpoint(integer,K,175)
    result={"status":"PASSED" if all(ep["checks"].values()) else "AUGMENTED ENDPOINT FAILED",
        "input_file":str(path),"input_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "flint_version":flint.__version__,"Arb_precision":ctx.prec,
        "local_primes_checked":count,"readback":saved,"rational_recursion":recur,
        "structural_checks":structure,"tilted_trace_Euler":KB,
        "endpoint_176":ep,"endpoint_175_diagnostic":ep_lower,
        "seconds_display_only":time.monotonic()-start}
    Path(args.output).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"Cjoint_integer":integer,
        "K_B_integer":K,"K_B_display_only":KB["K_B_display_only"],
        "176":ep["checks"],"176_ratio_display_only":ep["main_ratio_display_only"],
        "175":ep_lower["checks"],"structural_checks":structure,
        "seconds_display_only":result["seconds_display_only"]},indent=2),flush=True)
    if not all(ep["checks"].values()):
        raise SystemExit(1)


if __name__=="__main__":
    main()
