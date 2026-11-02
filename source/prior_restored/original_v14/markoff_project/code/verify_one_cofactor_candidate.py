#!/usr/bin/env python3
"""Read-only independent checks of a cofactor candidate and a frozen CRT run.

This verifies exact integer identities with FLINT arithmetic, not a new
determinant algorithm. The modular computation is a separate determinant
algorithm, but its unfinished CRT bound is not represented as a certificate.
"""
import argparse, datetime, functools, hashlib, json, math, pathlib, sys, time
sys.set_int_max_str_digits(0)
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'agents/computational_mechanisms/deps'))
from flint import fmpz

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def small_primes(n):
    a=bytearray(b'\x01')*(n+1); a[:2]=b'\x00\x00'
    for p in range(2,math.isqrt(n)+1):
        if a[p]: a[p*p::p]=b'\x00'*(((n-p*p)//p)+1)
    return [p for p in range(2,n+1) if a[p]]
def strip(n,ps):
    n=abs(n)
    if n==0:return 0
    for p in ps:
        while n%p==0:n//=p
    return n

ap=argparse.ArgumentParser()
ap.add_argument('--candidate',type=pathlib.Path,required=True)
ap.add_argument('--source',type=pathlib.Path,required=True)
ap.add_argument('--base',type=pathlib.Path,required=True)
ap.add_argument('--solve',type=pathlib.Path,required=True)
ap.add_argument('--crt',type=pathlib.Path,required=True)
ap.add_argument('--output',type=pathlib.Path,required=True)
args=ap.parse_args(); start=time.monotonic()
paths={k:getattr(args,k) for k in ['candidate','source','base','solve','crt']}
hashes={k:sha(v) for k,v in paths.items()}
c,s,b,u,crt=(read(paths[k]) for k in ['candidate','source','base','solve','crt'])
rows=s['rows']; n=len(rows[0]); i=c['binding']['column']; D=int(u['denominator'])
Q=list(map(int,u['numerators'])); determinant=int(b['minor_certificates'][0]['det'])
assert len(Q)==n and len(rows)>n and all(len(r)==n for r in rows)
assert c['binding']==u['binding'] and c['binding']['method_version']==2
assert c['binding']['equation']=='A Q = D e_i'
assert c['binding']['source_sha256']==hashes['source']==b['matrix_sha256']
assert c['binding']['base_minors_sha256_at_start']==hashes['base']
assert c['cofactor_checkpoint_sha256']==hashes['solve']
assert c['binding']['base_determinant_sha256']==hashlib.sha256(str(determinant).encode()).hexdigest()
assert b['minor_certificates'][0]['row_indices']==list(range(n))
assert int(c['base_determinant'])==determinant and int(c['denominator'])==D
assert D>0 and determinant%D==0 and functools.reduce(math.gcd,Q,D)==1
fq=list(map(fmpz,Q))
for j,row in enumerate(rows[:n]):
    total=sum((fmpz(a)*v for a,v in zip(row,fq)),fmpz(0))
    assert total==D*int(j==i),j
    if (j+1)%50==0: print(json.dumps({'exact_integer_rows_checked':j+1,'seconds':time.monotonic()-start}),flush=True)
ps=small_primes(2*s['nd']); roughs=[strip(determinant,ps)]
values=[]
for entry in c['replacement_minor_certificates']:
    ix=entry['row_indices']; assert len(ix)==n and len(set(ix))==n
    assert all(ix[j]==j for j in range(n) if j!=i)
    assert n<=ix[i]<len(rows)
    dot=int(sum((fmpz(a)*v for a,v in zip(rows[ix[i]],fq)),fmpz(0)))
    val=(determinant//D)*dot
    assert val==int(entry['det'])
    rough=strip(val,ps);assert rough==int(entry['rough_det'])
    roughs.append(rough);values.append((ix,val))
assert functools.reduce(math.gcd,roughs)==1==int(c['rough_gcd'])
assert crt['binding']['source_sha256']==hashes['source']
assert crt['binding']['shape']==[len(rows),n]
value=next(v for ix,v in values if ix==crt['binding']['row_indices'])
factor=int(crt['content_factor']); assert factor>0 and value%factor==0
reduced=value//factor
all_primes=set(); product=1
for rec in crt['records']:
    p,k,r=rec['p'],rec['k'],rec['r']
    assert k%2==1 and 0<k<2**r and p==k*2**r+1
    assert pow(rec['witness'],(p-1)//2,p)==p-1
    assert p not in all_primes;all_primes.add(p)
    assert reduced%p==rec['det_mod_p']
    product*=p
assert product==int(crt['modulus'])
assert reduced%product==int(crt['residue'])
assert all(sha(v)==hashes[k] for k,v in paths.items())
out={'status':'passed_exact_integer_identity_and_independent_modular_crosscheck',
     'd':s['d'],'nd':s['nd'],'shape':[len(rows),n],
     'inputs':{k:{'path':str(v),'sha256':hashes[k]} for k,v in paths.items()},
     'exact_equation':'A Q = D e_i','exact_rows_checked':n,
     'primitive_QD':True,'denominator_divides_base_determinant':True,
     'replacement_determinants_verified':len(values),'rough_gcd':'1',
     'independent_proth_moduli':len(all_primes),'crt_modulus_bits':product.bit_length(),
     'crt_hadamard_bound_bits':crt['hadamard']['bound_bits'],
     'crt_still_partial_not_an_exact_determinant_certificate':crt['status']!='complete',
     'arithmetic_engine':'FLINT exact integer multiplication/addition; separate modular determinant records',
     'base_determinant_imported_not_recomputed':True,
     'elapsed_seconds':time.monotonic()-start,
     'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out),flush=True)
