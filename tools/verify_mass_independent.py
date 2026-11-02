#!/usr/bin/env python3
"""Independent exact verifier of new mass weights, induction, and endpoint.

No producer module is imported. Rational expressions are rebuilt directly;
saved fixed-point quantities are checked as upper bounds, not trusted floats.
Inherited spectral nodes and curvature arguments retain their separate audit.
"""
from fractions import Fraction as F
from pathlib import Path
import os
import gzip, hashlib, json, sys, time
sys.set_int_max_str_digits(0)
ROOT=Path(os.environ['MARKOFF_WORKDIR'])
HERE=ROOT/'agents/analytic_mass_round3'
DATA=ROOT/'prior_restored/latest/positive157/complete'
start=time.monotonic()
def load(p):return json.loads(Path(p).read_text())
def f(x):return F(int(x['numerator']),int(x['denominator']))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def same_dependencies(c):
    for r in c['dependencies']:
        p=(ROOT/r['path'].split('/baragar_20261101_090837/',1)[1]); assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256']

mp=HERE/'MASS_MOMENT_CERTIFICATE.json'; ep=HERE/'MASS_ENDPOINT_CERTIFICATE.json'
c=load(mp); h=load(DATA/'parameters.json'); old=load(DATA/'positive157.json')
same_dependencies(c)
assert c['fixed_bits']==160 and c['integer_constant']==28924
S=1<<160
params=h['parameters']; meta=h['setup']; density=meta['density']
rows=params['active_primes']; assert len(rows)==74
assert [x['q'] for x in c['new_windows']]==[x['q'] for x in rows]
assert meta['column_prime_indices']==list(range(1,75))
pre_raw=gzip.decompress((DATA/'mellin.prefixes.json.gz').read_bytes())
node_meta=load(DATA/'mellin.json')
assert hashlib.sha256(pre_raw).hexdigest()==node_meta['prefix_raw_sha256']
prefix={r['panels']:r for r in json.loads(pre_raw)}
piL,piU=f(params['pi_lower']),f(params['pi_upper'])
slope=F(18,17); windows={}; fresh_windows=0
for i,(r,w) in enumerate(zip(rows,c['new_windows'])):
    q=r['q']; cap=f(meta['Euler_tail_caps'][i])
    if w['method']=='Euler':
        J=cap
    else:
        fresh_windows+=1
        lam=F(w['lambda']); k=w['split_panels']; N=w['end_panels']
        assert lam*density==N and 1<=k<=N
        lo=f(r['r_lower'])/q; hi=f(r['r_upper'])/q
        A=(1-lo)/slope; B=1+hi
        cc=piU*slope/(6*lam)
        phi=min(F(1),(1+hi)/2+(1-hi)*cc) if cc<F(1,2) else F(1)
        assert f(w['rho_lower'])==lo and f(w['rho_upper'])==hi
        for key,val in [('A_upper',A),('B_upper',B),('c_upper',cc),('phi_upper',phi)]:
            assert f(w[key])==val
        Sk=prefix[k]['pair_sums'][i]; Sn=prefix[N]['pair_sums'][i]
        harmonic=prefix[N]['harmonic_upper_scaled'][i]-prefix[k]['harmonic_lower_scaled'][i] if k<N else 0
        assert harmonic>=0
        val=(A*F(Sk,density*(1<<25))/piL+B*F(harmonic,1<<105)/piL
             +2*phi*F(Sn,density*(1<<25))/lam)*f(meta['curvature'][i]['Q_mellin_upper'])
        assert val==f(w['unrounded_window_upper'])
        J=min(cap,val)
    assert J==f(w['J_single_upper'])
    Jup=F(int(w['J_single_scaled']),S); assert J<=Jup<J+F(1,S)
    windows[q]=Jup

old_rows={r['q']:r for r in old['recursion_rows']}
new_rows=c['recursion_rows'];assert [r['q'] for r in new_rows]==[r['q'] for r in reversed(rows)]
C=F(1); tested=0; tails=0
for i,(r,cert) in enumerate(zip(reversed(rows),new_rows)):
    q=r['q']; ind=len(rows)-1-i; assert F(int(cert['C_next_scaled']),S)==C
    ru,wu,su=[f(r[k+'_upper']) for k in ['r','w','s']]
    cap=f(meta['Euler_tail_caps'][ind]); fm=f(r['Fmax_upper'])
    A=cert['first_omitted_a'];assert A>=1
    # This one inequality covers every omitted a>=A: 0<w<1,
    # B_q(a)^b<=Fmax and J<=cap, r^a<=1.
    assert 0<wu<1 and wu**A*fm*(cap+C)<=C;tails+=1
    assert [v['a'] for v in cert['candidates']]==list(range(1,A))
    by_a={v['a']:f(v['J_upper']) for v in old_rows[q]['candidates']}
    verified=[C]
    for v in cert['candidates']:
        a=v['a']; B=1+F(q-1,q)*sum((su**j for j in range(1,a+1)),F(0))
        Bup=F(int(v['B_scaled']),S); assert Bup>=B
        factor=F(int(v['factor_scaled']),S); assert factor**15>=Bup**2
        wp=F(int(v['w_power_scaled']),S);rp=F(int(v['r_power_scaled']),S)
        assert wp>=wu**a and rp>=ru**a
        geo=F(int(v['geometric_layer_sum_scaled']),S)
        assert geo>=sum((ru**j for j in range(a)),F(0))
        jm=F(int(v['J_mass_scaled']),S); jc=F(int(v['J_count_scaled']),S)
        assert jm>=geo*windows[q] and jc>=by_a.get(a,cap)
        J=F(int(v['J_selected_scaled']),S)
        # At least one independently valid branch is no larger than J.
        assert J>=min(jm,jc,cap)
        val=F(int(v['value_scaled']),S)
        assert val>=wp*factor*(J+rp*C)
        verified.append(val);tested+=1
    raw=F(int(cert['C_uncapped_scaled']),S); assert raw==max(verified)
    countcap=F(int(cert['C_count_cap_scaled']),S)
    assert countcap>=f(old_rows[q]['C_upper'])
    nxt=F(int(cert['C_mass_scaled']),S);assert nxt==min(raw,countcap) and nxt>=C
    C=nxt
assert C==f(c['exact_upper']) and C<28924

out={'status':'PASS_MASS_RECURSION','mass_constant':28924,'fresh_windows':fresh_windows,'finite_exponents':tested,'infinite_exponent_tails':tails,'script_sha256':sha(__file__)}
(ROOT/'root_logs/MASS_REPLAY_PUBLIC.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
