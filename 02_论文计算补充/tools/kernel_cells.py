#!/usr/bin/env python3
"""Rigorous CL majorant transport, fixed bandwidth Omega=32.

Uses the certified 15104-cell full-distribution Fourier table on [0,32].
This script does not regenerate that table or claim a new Markoff endpoint.
All new arithmetic uses Arb interval bounds; grid cells cover the continuum.
"""
from pathlib import Path
from fractions import Fraction as F
import json, hashlib, time, sys, platform
from flint import arb, acb, ctx

ctx.prec=192
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
TABLE=ROOT/'certificates/FREQUENCY_ENVELOPE.json'
TABLE_SHA=hashlib.sha256(TABLE.read_bytes()).hexdigest()
ETA=F(1,18)
OMEGA=32

def aa(q):
    q=F(q);return arb(q.numerator)/q.denominator
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def grid():
    return ([(F(k,16384),F(k+1,16384)) for k in range(4096)]
    +[(F(k,1024),F(k+1,1024)) for k in range(256,8192)]
    +[(F(k,128),F(k+1,128)) for k in range(1024,65536)])

def mhat_bound(l,u):
    delta=arb(OMEGA)/(2*arb.pi());lam=aa(ETA)/delta;r=(-lam).exp()
    ww=arb(aa((l+u)/2),aa((u-l)/2));v=ww/OMEGA
    z=acb(0,-2*arb.pi()*v).exp()
    num=1-v-acb(0,1)*lam*r*(1-z)/(2*arb.pi()*(1-r))
    den=delta*(1-r*z)
    return abs(num/den).upper()

def main():
    assert sha(TABLE)==TABLE_SHA
    table=json.loads(TABLE.read_text());cells=[]
    for i,(l,u) in enumerate(grid()):
        if l>=OMEGA:break
        assert u<=OMEGA
        weight=mhat_bound(l,u)*aa(u-l)/arb.pi()
        cells.append({'index':i,'l':str(l),'u':str(u),'fourier_weight_upper':str(weight.upper().fmpq())})
    assert len(cells)==15104 and F(cells[-1]['u'])==OMEGA
    out={'omega':OMEGA,'table_sha256':TABLE_SHA,'cells':cells}
    p=HERE/'bandlimited_kernel_cells.json';p.write_text(json.dumps(out,indent=2))
    (HERE/'BANDWIDTH_KERNEL_INPUTS.json').write_text(json.dumps({'kernel_cells_sha256':sha(p),'source_sha256':sha(Path(__file__))},indent=2))
    print('PASS_KERNEL_CELLS',len(cells))
if __name__=='__main__':main()
