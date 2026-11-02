#!/usr/bin/env python3
"""Verify a top-half Vandermonde determinant formula by exact arithmetic.

This diagnostic is conditional on the separately proved top-half formula;
optional --direct compares it with an ordinary integer determinant.
"""
import argparse,fractions,json,math,pathlib,time
from martin_kappa0 import ROOT,cheb_even,primes,smooth_strip
from flint import fmpz_poly,fmpz_mat


def formula(data):
    n=data['nd'];par=next(p for p in data['parameters'] if p['n']==n)
    G=fmpz_poly(par['g_T']);allroots=fmpz_poly(cheb_even(n-1));B=allroots//allroots.gcd(G)
    r=B.degree();assert r==par['m_count'];assert r-1<n/2
    disc=int(B.discriminant());res=int(B.resultant(G));B4=int(B(4))
    product_c=math.prod((-1)**j*math.comb(2*n-j-1,j) for j in range(r))
    D=fractions.Fraction(2*(-1)**n,n)**r*product_c*(-4)**(r*(r-1)//2)*disc*res*((-1)**r*B4)**(n-r+1)
    # Actual stored rows are m=r-1,...,0, while the formula uses m=0,...,r-1.
    D*=(-1)**(r*(r-1)//2)
    div=math.prod(int(meta['divided_small_content']) for meta in data['row_metadata'] if meta['n']==n)
    D/=div;assert D.denominator==1
    return int(D),{'r':r,'disc_B':str(disc),'resultant_B_G':str(res),'B_at_4':str(B4),'row_content_product':str(div)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('d',type=int,nargs='+');ap.add_argument('--direct',action='store_true');args=ap.parse_args();outpath=ROOT/'top_block_formula_results.json';out=json.loads(outpath.read_text()) if outpath.exists() else []
    for d in args.d:
        choices=list(ROOT.glob(f'stream_d{d}_matrix.json')) or list(ROOT.glob(f'*_d{d}_matrix.json'));src=choices[0];data=json.loads(src.read_text());st=time.time();D,parts=formula(data)
        result={'d':d,'nd':data['nd'],'source':src.name,'formula_seconds':time.time()-st,'formula_determinant':str(D),'determinant_rough':str(smooth_strip(D,primes(2*data['nd']))),**parts}
        if args.direct:
            r=parts['r'];A=fmpz_mat([row[-r:][::-1] for row in data['rows'][:r]]);t=time.time();direct=int(A.det());assert direct==D;result['direct_equal']=True;result['direct_seconds']=time.time()-t
        out=[old for old in out if old['d']!=d]+[result];out.sort(key=lambda z:z['d']);outpath.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['formula_determinant','disc_B','resultant_B_G','row_content_product']}),flush=True)


if __name__=='__main__':main()
