#!/usr/bin/env python3
"""Generic codimension-one certificate using n-1 shared source rows.

This applies to prime powers as well as primes.  It is elementary lattice
linear algebra, not a new Markoff-specific theorem.  It uses one (n-1)-order
determinant and one RHS instead of two n-order determinants.  If the common
row lattice is not smooth-saturated, this particular shortcut is
inconclusive even if the full matrix is good.
"""
import argparse,functools,hashlib,json,math,pathlib,time
from martin_kappa0 import ROOT,primes,smooth_strip
from flint import fmpz_mat,fmpq_mat,nmod_mat


def main():
    ap=argparse.ArgumentParser();ap.add_argument('d',type=int,nargs='+');args=ap.parse_args()
    for d in args.d:
        choices=list(ROOT.glob(f'stream_d{d}_matrix.json')) or list(ROOT.glob(f'*_d{d}_matrix.json'));src=choices[0];source=json.loads(src.read_text());rows=source['rows'];n=len(rows[0]);base_ids=list(range(1,n));base=[rows[i] for i in base_ids];st=time.time()
        R,rank=nmod_mat(base,1000000007).rref();assert rank==n-1
        pivots=[next(j for j in range(n) if int(R[i,j])) for i in range(n-1)];free=next(j for j in range(n) if j not in pivots);permutation=pivots+[free]
        A=fmpz_mat([[row[j] for j in pivots] for row in base]);b=fmpz_mat([[-row[free]] for row in base]);v=A.solve(b);assert fmpq_mat(A)*v==fmpq_mat(b)
        L=1
        for i in range(n-1):L=math.lcm(L,int(v[i,0].denominator))
        wp=[int(v[i,0].numerator)*(L//int(v[i,0].denominator)) for i in range(n-1)]+[L];w=[0]*n
        for j,k in enumerate(permutation):w[k]=wp[j]
        assert functools.reduce(math.gcd,w,0)==1 and all(sum(x*y for x,y in zip(row,w))==0 for row in base)
        D=int(A.det());assert D%L==0;g=abs(D)//L;ps=primes(2*source['nd']);rough=smooth_strip(g,ps)
        rec={'d':d,'nd':source['nd'],'source':src.name,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'shape':[len(rows),n],'method':'one determinant and a primitive cofactor vector for n-1 shared rows','base_row_indices':base_ids,'column_permutation':permutation,'primitive_kernel':list(map(str,w)),'kernel_denominator_lcm':str(L),'square_minor_determinant':str(D),'base_gcd_all_maximal_minors':str(g),'base_rough_gcd':str(rough),'base_seconds':time.time()-st}
        if rough!=1:rec['status']='shared_base_not_smooth';(ROOT/f'shared_cofactor_d{d}_certificate.json').write_text(json.dumps(rec,indent=2));print(json.dumps({'d':d,'status':rec['status'],'base_seconds':rec['base_seconds']}),flush=True);break
        residual=[];h=0
        for i,row in enumerate(rows):
            if i in base_ids:continue
            value=sum(x*y for x,y in zip(row,w));rv=smooth_strip(value,ps);residual.append({'row_index':i,'pairing':str(value),'rough_pairing':str(rv)});h=math.gcd(h,rv)
            if h==1:break
        rec.update({'residual_pairings':residual,'rough_gcd':str(h),'status':'certified_no_exception_p_above_bound' if h==1 else 'remaining_residual_gcd','total_seconds':time.time()-st})
        (ROOT/f'shared_cofactor_d{d}_certificate.json').write_text(json.dumps(rec,indent=2))
        print(json.dumps({k:v for k,v in rec.items() if k not in ['base_row_indices','column_permutation','primitive_kernel','kernel_denominator_lcm','square_minor_determinant','base_gcd_all_maximal_minors','residual_pairings']}),flush=True)


if __name__=='__main__':main()
