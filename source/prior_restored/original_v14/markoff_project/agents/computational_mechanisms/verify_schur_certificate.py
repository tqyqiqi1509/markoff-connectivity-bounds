#!/usr/bin/env python3
"""Independent builtins-only exact verification of Schur rank certificates."""
import fractions,hashlib,json,math,pathlib,sys,time
from verify_minor_certificates import bareiss,primes,strip


def verify(path):
    st=time.time();c=json.loads(path.read_text());src=path.parent/c['source'];data=src.read_bytes();assert hashlib.sha256(data).hexdigest()==c['source_sha256'];source=json.loads(data)
    ep=path.parent/c['elimination_file'];raw=ep.read_bytes();assert hashlib.sha256(raw).hexdigest()==c['elimination_sha256'];e=json.loads(raw)
    assert e['source_sha256']==c['source_sha256'] and e['columns_reversed']
    rows=[v[::-1] for v in source['rows']];r=e['top_block_size'];A=[v[:r] for v in rows[:r]];B=[v[r:] for v in rows[:r]];C=[v[:r] for v in rows[r:]];D=[v[r:] for v in rows[r:]];ps=primes(c['bound'])
    det=bareiss(A);assert str(det)==c['top_block_det'];assert strip(det,ps)==1
    X=[[fractions.Fraction(v) for v in row] for row in e['A_inverse_B']];k=len(X[0]);assert len(X)==r
    for i in range(r):
        for j in range(k):assert sum(A[i][t]*X[t][j] for t in range(r))==B[i][j]
    cols=list(map(int,e.get('column_scales',[]))) or [1]*k
    assert len(cols)==k and all(strip(q,ps)==1 for q in cols)
    target=e['normalized_schur_rows'];assert len(target)==len(C)
    for i in range(len(C)):
        scale=e['row_scales'][i];den=int(scale['denominator_lcm']);div=int(scale['divided_small_content']);assert strip(den,ps)==strip(div,ps)==1
        for j in range(k):
            S=(D[i][j]-sum(C[i][t]*X[t][j] for t in range(r)))*cols[j]
            assert S*den==target[i][j]*div
    g=0
    for minor in c['minor_certificates']:
        ids=minor['row_indices'];assert len(set(ids))==len(ids)==k
        v=bareiss([target[i] for i in ids]);assert str(v)==minor['det'];rough=strip(v,ps);assert str(rough)==minor['rough_det'];g=math.gcd(g,rough)
    assert g==1 and c['rough_gcd']=='1'
    return {'certificate':path.name,'verified':True,'method':'independent Python Fraction matrix identities and Bareiss determinants','d':c['d'],'top_block':r,'schur_columns':k,'seconds':time.time()-st}


if __name__=='__main__':
    for arg in sys.argv[1:]:print(json.dumps(verify(pathlib.Path(arg))),flush=True)
