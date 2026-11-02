#!/usr/bin/env python3
"""Exhaustive exact-rational structure/hash check, not a numerical replay.

No imports from the producer; all interval coverage and saved upper-bound
composition are checked using Python Fraction only. The root separately
recomputes transcendental arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read(name): return json.loads((HERE/name).read_text())
def upper(pack):
    lo,hi=F(pack['lower']),F(pack['upper']); assert lo<=hi
    return hi

table=read('EULER_FREQUENCY_CERTIFICATE.json')
cert=read('FIXED_GAUSSIAN_TRACE_CERTIFICATE.json')
gate=read('FIXED_PARAMETER_GATES.json')
assert table['parameters']==cert['parameters']==gate['parameters']
assert table['producer_sha256']==cert['producer_sha256']==gate['producer_sha256']==sha(HERE/'certify_fixed_gaussian.py')
assert cert['euler_table_sha256']==sha(HERE/'EULER_FREQUENCY_CERTIFICATE.json')
assert cert['fixed_gates_sha256']==table['fixed_gates_sha256']==sha(HERE/'FIXED_PARAMETER_GATES.json')
assert table['status']=='PASS_ALL_PRIMES_EXPONENTS_CONTINUOUS_FREQUENCY_BOXES'
assert cert['status']=='PASS_FIXED_GAUSSIAN_TRACE_ALL_PRIMES_EXPONENTS_FREQUENCIES_VARIANCES_PENDING_ROOT_REVIEW'
for entry in cert['dependencies']:
    assert sha(ROOT/entry['path'])==entry['sha256']
    assert (ROOT/entry['path']).stat().st_size==entry['bytes']
assert F(gate['D_expression']['lower'])>F(1693,200)
lattice=gate['lattice_variance_gate']
assert F(lattice['log2']['lower'])>F(693,1000)
coef=F(43,1000)-F(61,10000)*F(693,2000)
assert F(lattice['coefficient_upper'])==coef>0
assert F(lattice['D_minus_coefficient_times_207'])==F(1693,200)-207*coef>0

sieve=bytearray(b'\x01')*4099; sieve[:2]=b'\x00\x00'
for n in range(2,math.isqrt(4098)+1):
    if sieve[n]:
        for j in range(n*n,4099,n):sieve[j]=0
primes=[q for q in range(2,4099) if sieve[q]]
assert [r['q'] for r in table['local_tails']]==primes
assert len(primes)==table['prime_count']==564
assert sum(r['first_omitted_exponent']-1 for r in table['local_tails'])==table['finite_nonzero_exponents']==1355
assert upper(table['prime_tail_guard'])<1
assert all(upper(r['tail_guard'])<1 for r in table['local_tails'])
grids=([(F(i,16384),F(i+1,16384)) for i in range(4096)]
       +[(F(i,1024),F(i+1,1024)) for i in range(256,8192)]
       +[(F(i,128),F(i+1,128)) for i in range(1024,65536)])
assert len(grids)==table['frequency_grid']['total_boxes']==len(table['frequency_upper_values'])==76544
assert grids[0][0]==0 and grids[-1][1]==512
assert all(a[1]==b[0] for a,b in zip(grids,grids[1:]))
assert all(0<F(s)<=1 for s in table['frequency_upper_values'])

edges=[F(i,16) for i in range(3312,4097)]+[F(i) for i in range(257,1601)]
rows=cert['variance_rows']; assert len(rows)==cert['variance_boxes']==2128
tail=upper(cert['frequency_infinite_tail_upper'])
highA=upper(cert['finite_high_frequency_A_upper'])
farG=upper(cert['finite_high_frequency_G_treatment']['far_uniform_G_upper'])
maxval=F(0);maxbox=None;usedBE=0
for row,vl,vu in zip(rows,edges,edges[1:]):
    assert list(map(F,row['V_interval']))==[vl,vu]
    tm=(F(43,1000)*vu-F(1693,200))/F(61,10000)
    assert F(row['T_upper']['lower'])<=tm<=F(row['T_upper']['upper'])
    for name in ['Gaussian_root_interval','BE_normal_root_interval']:
        l,h=map(F,row[name]);assert l<h
    finite=upper(row['finite_error_upper'])
    assert finite>=upper(row['finite_low_error_upper'])+highA+upper(row['finite_near_high_G_upper'])+farG
    rep=upper(row['replacement_upper'])
    assert rep>=upper(row['Gaussian_peak_upper'])+finite+tail
    arg=min(vu,max(vl,3*F(1693,200)/F(43,1000)))
    assert F(row['BE_error_argmax'])==arg
    be=upper(row['BE_upper'])
    assert be>=upper(row['BE_normal_H_upper'])+upper(row['BE_error_upper'])
    full=upper(row['total_upper'])
    assert full>=min(rep,be)
    assert full<F(189,500)
    if row['chosen_complete_bound']=='BE':
        usedBE+=1;assert be<rep
    else:assert row['chosen_complete_bound']=='Gaussian_replacement' and rep<=be
    if full>maxval:maxval=full;maxbox=row['V_interval']
assert maxval==upper(cert['maximum_upper']) and maxbox==cert['maximum_box']
assert usedBE==cert['boxes_using_BE']
infinite=upper(cert['infinite_variance_tail_upper'])
assert infinite<F(189,500) and cert['infinite_variance_tail_from']==1600
assert upper(cert['global_upper'])==max(maxval,infinite)
assert F(cert['gamma_margin']['lower'])>0
out={'status':'PASS_EXHAUSTIVE_EXACT_RATIONAL_STRUCTURE_AND_BINDINGS',
     'frequency_values':76544,'continuous_variance_boxes':2128,
     'all_prime_tail_rows':564,'all_nonzero_candidate_exponents':1355,
     'maximum_box':maxbox,'maximum_saved_upper':str(maxval),
     'global_saved_upper':str(upper(cert['global_upper'])),'boxes_using_BE':usedBE,
     'producer_sha256':sha(HERE/'certify_fixed_gaussian.py'),
     'frequency_sha256':sha(HERE/'EULER_FREQUENCY_CERTIFICATE.json'),
     'variance_sha256':sha(HERE/'FIXED_GAUSSIAN_TRACE_CERTIFICATE.json'),
     'checker_sha256':sha(Path(__file__)),
     'scope':'All saved rational values, coverage and compositions checked. This does not independently verify transcendental enclosures; independent arithmetic replay belongs to root.'}
(HERE/'EXACT_STRUCTURE_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
