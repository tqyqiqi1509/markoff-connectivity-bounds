#!/usr/bin/env python3
"""Four-block rank certificate through a codimension-one lower block.

For odd prime d>=5:
  lower rows n=d,2d,3d: (3d-3) x (3d-2), saturated away from <=6d;
  fourth block n=4d: r=2d-2 rows, with a smooth invertible top-half block A.
The sole lower-block quotient functional is its primitive kernel w.
Solve A*z = fourth_block_low_columns*w.  Full column rank modulo p>8d is
equivalent to some z_i, i>=d, being nonzero.  Therefore a scalar gcd of
these numerators completes the certificate.
"""
import argparse,functools,hashlib,json,math,pathlib,time
from martin_kappa0 import ROOT,factor,primes,smooth_strip
from top_block_determinant_formula import formula
from flint import fmpz_mat,fmpq_mat


def lower_kernel(source,src):
    d=source['d'];selected=[(row,meta) for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']<=3*d]
    restored=[smooth_strip(int(meta['divided_small_content']),primes(6*d)) for row,meta in selected]
    rows=[[x*f for x in row[:3*d-2]] for (row,meta),f in zip(selected,restored)]
    sha=hashlib.sha256(src.read_bytes()).hexdigest();path=ROOT/f'three_block_d{d}_kernel.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('source_sha256')==sha:
            w=list(map(int,old['primitive_kernel']));assert all(sum(x*y for x,y in zip(v,w))==0 for v in rows)
            return old,rows
    st=time.time();A=fmpz_mat([row[:-1] for row in rows]);b=fmpz_mat([[-row[-1]] for row in rows]);t=time.time();x=A.solve(b);solve_seconds=time.time()-t;assert fmpq_mat(A)*x==fmpq_mat(b)
    L=1
    for i in range(x.nrows()):L=math.lcm(L,int(x[i,0].denominator))
    w=[int(x[i,0].numerator)*(L//int(x[i,0].denominator)) for i in range(x.nrows())]+[L]
    assert functools.reduce(math.gcd,w,0)==1
    t=time.time();det=int(A.det());det_seconds=time.time()-t;assert det%L==0;g=abs(det)//L
    obj={'d':d,'N':3*d,'source':src.name,'source_sha256':sha,'shape':[len(rows),len(rows[0])],'restored_content_factors':list(map(str,restored)),'primitive_kernel':list(map(str,w)),'last_coordinate':str(L),'first_square_minor_det':str(det),'gcd_all_maximal_minors':str(g),'gcd_rough_above_6d':str(smooth_strip(g,primes(6*d))),'solve_seconds':solve_seconds,'determinant_seconds':det_seconds,'seconds':time.time()-st}
    path.write_text(json.dumps(obj,indent=2));return obj,rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('d',type=int,nargs='+');ap.add_argument('--direct-top-det',action='store_true');args=ap.parse_args()
    for d in args.d:
        assert d>=5 and factor(d)==[(d,1)]
        src=ROOT/f'stream_d{d}_matrix.json'
        if not src.exists():src=next(ROOT.glob(f'*_d{d}_matrix.json'))
        st=time.time();source=json.loads(src.read_text());assert source['nd']==4*d
        lower,_=lower_kernel(source,src);w=list(map(int,lower['primitive_kernel']))
        rec={'d':d,'nd':4*d,'bound':8*d,'source':src.name,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'method':'three-block saturated codimension-one kernel plus top-half scalar residual gcd','lower_certificate':lower,'lower_seconds_this_run':time.time()-st}
        print(json.dumps({k:v for k,v in rec.items() if k!='lower_certificate'}),flush=True)
        if lower['gcd_rough_above_6d']!='1':
            rec['status']='lower_block_obstruction';(ROOT/f'codim_d{d}_certificate.json').write_text(json.dumps(rec,indent=2));break
        top=[row for row,meta in zip(source['rows'],source['row_metadata']) if meta['n']==4*d];r=len(top);assert r==2*d-2
        A=fmpz_mat([row[-r:][::-1] for row in top]);rhs=[sum(x*y for x,y in zip(row[:3*d-2],w)) for row in top]
        t=time.time();z=A.solve(fmpz_mat([[v] for v in rhs]));assert fmpq_mat(A)*z==fmpq_mat([[v] for v in rhs]);rec['top_single_rhs_solve_seconds']=time.time()-t
        topdet,parts=formula(source);assert smooth_strip(topdet,primes(8*d))==1
        rec.update({'top_block_size':r,'top_block_determinant':str(topdet),'top_block_formula_parts':parts,'top_rhs':list(map(str,rhs)),'top_solution':[str(z[i,0]) for i in range(r)],'residual_indices':list(range(d,r))})
        if args.direct_top_det:
            t=time.time();assert int(A.det())==topdet;rec['direct_top_det_verified']=True;rec['direct_top_det_seconds']=time.time()-t
        ps=primes(8*d);g=0
        for i in range(r):assert smooth_strip(int(z[i,0].denominator),ps)==1
        for i in range(d,r):g=math.gcd(g,smooth_strip(int(z[i,0].numerator),ps))
        rec['rough_gcd']=str(g);rec['status']='certified_no_exception_p_above_bound' if g==1 else 'residual_rough_gcd';rec['total_seconds']=time.time()-st
        path=ROOT/f'codim_d{d}_certificate.json';path.write_text(json.dumps(rec,indent=2))
        print(json.dumps({k:v for k,v in rec.items() if k not in ['lower_certificate','top_block_determinant','top_block_formula_parts','top_rhs','top_solution','residual_indices']}),flush=True)
        if g!=1:
            (ROOT/'global_certificate_stop_requested').write_text(f'Codimension-one full-rank certificate has a residual rough gcd for d={d}; inspect {path.name} before any further expansion.\n')
            break


if __name__=='__main__':main()
