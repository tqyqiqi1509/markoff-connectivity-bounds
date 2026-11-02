#!/usr/bin/env python3
"""Exact two-dimensional Newton recurrence for Phi at kappa=0."""
import argparse, functools, json, math, pathlib, resource, time
from martin_kappa0 import fmpz_poly, ReductionCache, independent_reduction

class NewtonReductionCache:
    def __init__(self,maxa,progress=True):
        self.pairs={};self.maxa=maxa;self.start=time.time()
        self.coefficients={b:[0]+[(-1)**k*(2*b*math.comb(2*b-k,k)//(2*b-k)) for k in range(1,b+1)] for b in range(1,maxa//2+1)}
        for s in range(maxa+1):
            for b in range(s//2+1):
                a=s-b
                if b==0:out=fmpz_poly([0]*a+[1])
                else:
                    out=2*self[a-b,b]
                    co=self.coefficients[b]
                    for k in range(1,b+1):
                        inner=fmpz_poly([])
                        for j in range(k+1):inner+=math.comb(k,j)*self[a-k+j,b-j]
                        out-=co[k]*inner
                self.pairs[a,b]=out
            if progress and (s%10==0 or s==maxa):print(json.dumps({'stage':'newton_cache','sum':s,'max':maxa,'elapsed':time.time()-self.start,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}),flush=True)
    def __getitem__(self,ab):
        a,b=ab
        return self.pairs[(max(a,b),min(a,b))]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--maxsum',type=int,default=30);ap.add_argument('--original-check',type=int,default=30);args=ap.parse_args()
    st=time.time();new=NewtonReductionCache(args.maxsum)
    tnew=time.time()-st
    st=time.time();old=ReductionCache(args.maxsum,False);told=time.time()-st
    checks=0
    for key,out in new.pairs.items():assert out==old[key],key;checks+=1
    original=0;st=time.time()
    for a,b in new.pairs:
        if a+b>args.original_check:continue
        q=independent_reduction(2*a,2*b,0)
        p=fmpz_poly([q[2*i] for i in range(q.degree()//2+1)])
        assert p==new[a,b],(a,b)
        original+=1
    rec={'kappa':0,'max_half_total_degree':args.maxsum,'pair_checks_against_even_three_dimensional':checks,'checks_against_original_one_step_recurrence':original,'original_check_max_half_total_degree':args.original_check,'newton_seconds':tnew,'even_3D_seconds':told,'original_check_seconds':time.time()-st,'all_equal':True,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(rec),flush=True)
    pathlib.Path(f'newton_validation_{args.maxsum}.json').write_text(json.dumps(rec,indent=2))

if __name__=='__main__':main()
