#!/usr/bin/env python3
"""Independent verification of the codimension-one rank certificate.

Default: Python integers/Fraction and independent Bareiss determinants.
For large certificates --flint changes only the determinant implementation.
No top-half formula is assumed by the verifier: it evaluates the top
integer determinant directly from the source matrix.
"""
import argparse,fractions,functools,hashlib,json,math,pathlib,sys,time
from verify_minor_certificates import bareiss,primes,strip


def verify(path,determinant):
    st=time.time();c=json.loads(path.read_text());src=path.parent/c['source'];raw=src.read_bytes();assert hashlib.sha256(raw).hexdigest()==c['source_sha256'];source=json.loads(raw)
    d=c['d'];assert d>=5 and all(d%q for q in range(2,math.isqrt(d)+1))
    assert source['d']==d and source['kappa']==0
    assert source['nd']==c['nd']==4*d;assert c['bound']==8*d
    assert len(source['rows'])==len(source['row_metadata'])==5*d-5
    assert all(len(row)==4*d-2 for row in source['rows'])
    # The quotient argument needs exact support, not merely the labels n.
    # The generator already enforces this degree bound before serializing.
    assert all(all(v==0 for v in row[meta['n']-2:])
               for row,meta in zip(source['rows'],source['row_metadata']))
    ps=primes(8*d)
    for meta in source['row_metadata']:assert strip(int(meta['divided_small_content']),ps)==1
    lower=c['lower_certificate'];assert lower['source_sha256']==c['source_sha256'];w=list(map(int,lower['primitive_kernel']));L=int(lower['last_coordinate']);assert w[-1]==L>0 and functools.reduce(math.gcd,w,0)==1
    selected=[(row,meta) for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']<=3*d];restore=[strip(int(meta['divided_small_content']),primes(6*d)) for row,meta in selected];assert list(map(str,restore))==lower['restored_content_factors']
    rows=[[x*f for x in row[:3*d-2]] for (row,meta),f in zip(selected,restore)]
    assert len(rows)==3*d-3 and len(w)==3*d-2
    assert all(sum(x*y for x,y in zip(row,w))==0 for row in rows)
    D=determinant([row[:-1] for row in rows]);assert str(D)==lower['first_square_minor_det'];assert D%L==0;g=abs(D)//L;assert str(g)==lower['gcd_all_maximal_minors'];assert strip(g,primes(6*d))==1
    top=[row for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']==4*d];r=len(top);assert r==2*d-2==c['top_block_size'];A=[row[-r:][::-1] for row in top]
    det=determinant(A);assert str(det)==c['top_block_determinant'];assert strip(det,ps)==1
    rhs=[sum(x*y for x,y in zip(row[:3*d-2],w)) for row in top];assert list(map(str,rhs))==c['top_rhs']
    z=[fractions.Fraction(v) for v in c['top_solution']];assert len(z)==r and all(strip(v.denominator,ps)==1 for v in z)
    assert all(sum(a*b for a,b in zip(row,z))==s for row,s in zip(A,rhs));assert c['residual_indices']==list(range(d,r))
    h=functools.reduce(math.gcd,(strip(v.numerator,ps) for v in z[d:]),0);assert h==1 and c['rough_gcd']=='1'
    return {'certificate':path.name,'verified':True,'d':d,'lower_shape':[3*d-3,3*d-2],'top_block_size':r,'scalar_residual_count':d-2,'seconds':time.time()-st}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('files',nargs='+');ap.add_argument('--flint',action='store_true');args=ap.parse_args();determinant=bareiss
    if args.flint:
        sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent/'deps'));from flint import fmpz_mat
        determinant=lambda A:int(fmpz_mat(A).det())
    for fn in args.files:
        result=verify(pathlib.Path(fn),determinant);result['method']='Python Fraction identities; '+('FLINT determinants' if args.flint else 'independent pure-Python Bareiss determinants');print(json.dumps(result),flush=True)
