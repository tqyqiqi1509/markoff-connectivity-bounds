#!/usr/bin/env python3
"""Verify Martin kappa=0 determinant certificates.

Default determinant routine is pure-Python fraction-free Bareiss, independent
of the generator's FLINT determinant. --flint is available for larger files.
This verifies finite arithmetic certificates, not the cited mathematical theorem.
"""
import argparse, hashlib, json, math, pathlib, sys, time
sys.set_int_max_str_digits(0)

def bareiss(A):
    A=[list(map(int,row)) for row in A];n=len(A);sgn=1;prev=1
    if n==0:return 1
    for k in range(n-1):
        if not A[k][k]:
            pivot=next((i for i in range(k+1,n) if A[i][k]),None)
            if pivot is None:return 0
            A[k],A[pivot]=A[pivot],A[k];sgn=-sgn
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                quotient,remainder=divmod(numerator,prev)
                if remainder:raise ArithmeticError('Bareiss non-exact division')
                A[i][j]=quotient
            A[i][k]=0
        prev=pivot
    return sgn*A[-1][-1]

def primes(n):return [q for q in range(2,n+1) if all(q%t for t in range(2,math.isqrt(q)+1))]
def strip(x,ps):
    x=abs(x)
    if not x:return x
    for p in ps:
        while x%p==0:x//=p
    return x

def main():
    ap=argparse.ArgumentParser();ap.add_argument('certificates',nargs='+');ap.add_argument('--flint',action='store_true');args=ap.parse_args()
    if args.flint:
        sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent/'deps'))
        from flint import fmpz_mat
        determinant=lambda A:int(fmpz_mat(A).det())
    else:determinant=bareiss
    for fn in args.certificates:
        ts=time.time();cf=pathlib.Path(fn);c=json.loads(cf.read_text());mf=cf.parent/c['matrix_file'];data=mf.read_bytes()
        assert hashlib.sha256(data).hexdigest()==c['matrix_sha256']
        obj=json.loads(data);rows=obj['rows'];n=obj['nd']-2;ps=primes(2*obj['nd'])
        assert all(len(row)==n for row in rows)
        for m in obj['row_metadata']:assert strip(int(m['divided_small_content']),ps)==1
        if 'codimension_one_certificate_file' in c:
            nested=cf.parent/c['codimension_one_certificate_file'];assert hashlib.sha256(nested.read_bytes()).hexdigest()==c['codimension_one_certificate_sha256']
            from verify_codimension_one_certificate import verify
            result=verify(nested,determinant);result.update({'wrapper_certificate':cf.name,'arithmetic_verified':True,'method':'codimension_one; '+('FLINT determinants' if args.flint else 'independent pure-Python Bareiss determinants')})
            print(json.dumps(result),flush=True);continue
        g=0
        for m in c['minor_certificates']:
            indices=m['row_indices'];assert len(set(indices))==len(indices)==n
            det=determinant([rows[i] for i in indices]);assert str(det)==m['det']
            rough=strip(det,ps);assert str(rough)==m['rough_det'];g=math.gcd(g,rough)
        assert str(g)==c['rough_gcd'];assert g==1
        print(json.dumps({'certificate':cf.name,'matrix_sha256':c['matrix_sha256'],'d':c['d'],'arithmetic_verified':True,'method':'FLINT' if args.flint else 'independent_pure_python_Bareiss','minor_count':len(c['minor_certificates']),'seconds':time.time()-ts}),flush=True)

if __name__=='__main__':main()
