#!/usr/bin/env python3
"""Exact full source-row comparison with the old forward polynomial engine.

The main engine propagates source coefficients backward through the DAG.
This verifier builds the forward polynomial E(a,b,c) values and combines
them using an independently assembled general-d source polynomial.
"""
import datetime,hashlib,json,pathlib,resource,sys,time
from martin_kappa0 import ROOT,matrix_for
from composite_parameters import parameters
from flint import fmpz_poly
sys.set_int_max_str_digits(0)

class ForwardPolynomialCache:
    """Same exact forward recurrence, with active-memory rather than cache guard."""
    def __init__(self,maxa):
        self.pairs={};previous={};start=time.time()
        for degree in range(maxa+1):
            current={}
            def E(a,b,c):return previous[tuple(sorted((a,b,c),reverse=True))]
            for b in range(degree//2+1):
                for c in range(min(b,degree-2*b)+1):
                    a=degree-b-c
                    if b==0:value=fmpz_poly([0]*a+[1])
                    elif c==0:value=2*(E(a,b-1,0)+E(a-1,b,0)+E(a-1,b-1,1))
                    else:value=E(a+1,b-1,c-1)+E(a-1,b+1,c-1)+E(a-1,b-1,c+1)+2*(E(a,b,c-1)+E(a,b-1,c)+E(a-1,b,c))
                    current[a,b,c]=value
                    if c==0:self.pairs[a,b]=value
            previous=current
            if degree%20==0 or degree==maxa:
                rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                stat={line.split()[0]:int(line.split()[1]) for line in pathlib.Path('/sys/fs/cgroup/memory.stat').read_text().splitlines()}
                print(json.dumps({'stage':'forward_polynomial_cache','sum':degree,'max':maxa,'elapsed':time.time()-start,'rss_kib':rss,'shared_anon_bytes':stat['anon'],'reclaimable_file_bytes':stat['file']}),flush=True)
                if rss>1000*1024 or stat['anon']>6*1024**3:raise MemoryError('active-memory diagnostic guard')
    def __getitem__(self,key):
        a,b=key;return self.pairs[max(a,b),min(a,b)]

d=15;src=ROOT/f'composite_d{d}_matrix.json';data=json.loads(src.read_text());nd=data['nd'];_,params=parameters(d,nd)
maxa=max(2*n+(len(g)-1)+m-1 for n,m,g in params);st=time.time()
cache=ForwardPolynomialCache(maxa)
rows,meta=matrix_for(d,nd,params,cache)
assert rows==data['rows'] and meta==data['row_metadata']
# The three exact low-coefficient constraints uniquely reconstruct a0,a1,a2
# from the retained coefficients; both engines explicitly checked them.
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'d':d,'nd':nd,'rows':len(rows),'retained_coefficient_checks':sum(len(r) for r in rows),'full_source_coefficients_including_low_terms':sum(m['n']+1 for m in meta),'all_retained_coefficients_and_normalizations_equal':True,'all_low_coefficients_follow_from_verified_exact_identities':True,'source':src.name,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'method':'old forward exact polynomial recurrence vs new backward scalar coefficient propagation; Python monic-lcm filter vs C++ inclusion-exclusion filter','max_half_degree':maxa,'seconds':time.time()-st,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(ROOT/'composite_d15_forward_polynomial_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
