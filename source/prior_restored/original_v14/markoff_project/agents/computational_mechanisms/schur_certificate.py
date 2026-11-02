#!/usr/bin/env python3
"""Peel the highest n block by an exact smooth Schur complement.

The selected top block A is checked to have determinant containing only
primes <= 2*n_d.  Its Schur complement therefore has exactly the remaining
rank modulo every larger prime.  Rational arithmetic is exact; each
denominator and each removed row-content divisor is separately checked to
be smooth.  The complete elimination matrix X=A^-1 B is saved for a verifier
to check A*X=B without trusting a linear-system solver.
"""
import argparse, functools, hashlib, json, math, pathlib, resource, sys, time
from martin_kappa0 import ROOT,primes,smooth_strip,gcd_minors
from flint import fmpz_mat,fmpq_mat


def matrix_rows(M):
    return [[str(M[i,j]) for j in range(M.ncols())] for i in range(M.nrows())]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('d',type=int);ap.add_argument('--source');ap.add_argument('--maxminors',type=int,default=32);ap.add_argument('--column-denominators',action='store_true');args=ap.parse_args()
    src=pathlib.Path(args.source) if args.source else next(ROOT.glob(f'stream_d{args.d}_matrix.json'))
    original=json.loads(src.read_text());rows=[r[::-1] for r in original['rows']];nd=original['nd'];bound=2*nd;ps=primes(bound)
    r=sum(meta['n']==nd for meta in original['row_metadata']);m=len(rows);n=len(rows[0]);st=time.time()
    assert all(meta['n']==nd for meta in original['row_metadata'][:r])
    A=fmpz_mat([v[:r] for v in rows[:r]]);B=fmpz_mat([v[r:] for v in rows[:r]])
    C=fmpz_mat([v[:r] for v in rows[r:]]);D=fmpz_mat([v[r:] for v in rows[r:]])
    t=time.time();det=int(A.det());rough=smooth_strip(det,ps)
    rec={'d':args.d,'nd':nd,'bound':bound,'source':src.name,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'original_shape':[m,n],'top_block_size':r,'top_block_det':str(det),'top_block_det_rough':str(rough),'top_block_det_seconds':time.time()-t}
    print(json.dumps({k:v for k,v in rec.items() if k!='top_block_det'}),flush=True)
    if rough!=1:
        rec['status']='top_block_not_smooth';(ROOT/f'schur_d{args.d}_certificate.json').write_text(json.dumps(rec,indent=2));return
    t=time.time();X=A.solve(B);assert fmpq_mat(A)*X==fmpq_mat(B);rec['solve_seconds']=time.time()-t
    t=time.time();S=fmpq_mat(D)-fmpq_mat(C)*X;rec['schur_multiply_seconds']=time.time()-t
    out=[];scales=[];column_scales=[]
    if args.column_denominators:
        for j in range(S.ncols()):
            den=1
            for i in range(S.nrows()):den=math.lcm(den,int(S[i,j].denominator))
            assert smooth_strip(den,ps)==1
            column_scales.append(str(den))
            for i in range(S.nrows()):S[i,j]*=den
    for i in range(S.nrows()):
        den=1
        for j in range(S.ncols()):den=math.lcm(den,int(S[i,j].denominator))
        assert smooth_strip(den,ps)==1
        row=[int(S[i,j].numerator)*(den//int(S[i,j].denominator)) for j in range(S.ncols())]
        co=functools.reduce(math.gcd,row,0);div=abs(co)//smooth_strip(co,ps) if co else 1
        out.append([x//div for x in row]);scales.append({'denominator_lcm':str(den),'divided_small_content':str(div)})
    prefix=f'schur_{"columns_" if args.column_denominators else ""}d{args.d}'
    elim={'source':src.name,'source_sha256':rec['source_sha256'],'columns_reversed':True,'top_block_size':r,'A_inverse_B':matrix_rows(X),'normalized_schur_rows':out,'row_scales':scales,'column_scales':column_scales}
    fp=ROOT/f'{prefix}_elimination.json';fp.write_text(json.dumps(elim,separators=(',',':')))
    rec.update({'elimination_file':fp.name,'elimination_sha256':hashlib.sha256(fp.read_bytes()).hexdigest(),'schur_shape':[len(out),len(out[0])],'schur_entry_max_bits':max(abs(x).bit_length() for row in out for x in row),'reduction_seconds':time.time()-st})
    print(json.dumps({k:v for k,v in rec.items() if k!='top_block_det'}),flush=True)
    rec.update(gcd_minors(out,bound,args.maxminors));rec['total_seconds']=time.time()-st;rec['maxrss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (ROOT/f'{prefix}_certificate.json').write_text(json.dumps(rec,indent=2))
    print(json.dumps({k:v for k,v in rec.items() if k not in ['top_block_det','minor_certificates']}),flush=True)


if __name__=='__main__':main()
