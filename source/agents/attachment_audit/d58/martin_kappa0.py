#!/usr/bin/env python3
"""Exact kappa=0 specialization of Martin Algorithm 1.

Arithmetic: python-flint fmpz_poly/fmpz_mat. No floating rank evidence is used.
This is a sufficiency certificate conditional on audited Martin Theorem 7.1,
not an independent proof of that theorem. Matrices retain columns x^6,...,x^(2n).
"""
import argparse, functools, hashlib, json, math, pathlib, resource, sys, time
sys.set_int_max_str_digits(0)
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'deps'))
from flint import fmpz_poly, fmpz_mat

def factor(n):
    fs=[];q=2
    while q*q<=n:
        if n%q==0:
            e=0
            while n%q==0:n//=q;e+=1
            fs.append((q,e))
        q+=1
    if n>1:fs.append((n,1))
    return fs

def divisors(n):
    ds=[1]
    for p,e in factor(n):ds=[d*p**i for d in ds for i in range(e+1)]
    return ds

def phi(n):
    r=n
    for p,e in factor(n):r=r//p*(p-1)
    return r

def val(n,p):
    e=0
    while n%p==0:n//=p;e+=1
    return e

def primes(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,math.isqrt(p)+1))]

def smooth_strip(n,ps):
    if not n:return 0
    n=abs(n)
    for p in ps:
        while n%p==0:n//=p
    return n

def cheb_even(k):
    """U_k(x/2), multiplied by x iff k odd; then written in T=x^2."""
    a=fmpz_poly([1]);b=fmpz_poly([0,1]);x=b
    if k==0:q=a
    elif k==1:q=b
    else:
        for j in range(2,k+1):a,b=b,x*b-a
        q=b
    if k%2:q*=x
    assert all(q[i]==0 for i in range(1,q.degree()+1,2))
    return [int(q[i]) for i in range(0,q.degree()+1,2)]

def parameters(d,nd=None):
    fs=factor(d)
    if len(fs)!=1 or d<4:raise ValueError('d must be prime power >=4')
    q,e=fs[0]
    if nd is None:nd=d*({2:6,3:5}.get(q,4))
    if nd<3*d:raise ValueError('Theorem 7.1 needs nd>=3d')
    params=[]
    for n in range(d,nd+1,d):
        raw=sum(phi(2*n//i) for i in divisors(2*n//(d*math.gcd(2,q))));mn=(raw+3)//4
        k=n*d//q**(val(n,q)+1)-1
        g=cheb_even(k)
        params.append((n,mn,g))
    return nd,params

class ReductionCache:
    def __init__(self,maxa,progress=True,rss_limit_mb=1400):
        self.pairs={};self.prev={};self.maxa=-1;self.rss_limit_mb=rss_limit_mb
        self.extend(maxa,progress)
    def extend(self,maxa,progress=True):
        prev=self.prev;self.prev={};start=time.time()
        for s in range(self.maxa+1,maxa+1):
            cur={}
            def get(a,b,c):return prev[tuple(sorted((a,b,c),reverse=True))]
            for b in range(s//2+1):
                for c in range(min(b,s-2*b)+1):
                    a=s-b-c
                    if b==0:
                        out=fmpz_poly([0]*a+[1])
                    elif c==0:
                        out=2*(get(a,b-1,0)+get(a-1,b,0)+get(a-1,b-1,1))
                    else:
                        out=(get(a+1,b-1,c-1)+get(a-1,b+1,c-1)+get(a-1,b-1,c+1)
                             +2*(get(a,b,c-1)+get(a,b-1,c)+get(a-1,b,c)))
                    cur[(a,b,c)]=out
                    if c==0:self.pairs[(a,b)]=out
            prev=cur
            rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            if progress and (s%20==0 or s==maxa):print(json.dumps({'stage':'reduction_cache','sum':s,'max':maxa,'states_layer':len(cur),'elapsed':round(time.time()-start,3),'maxrss_kib':rss}),flush=True)
            if rss>self.rss_limit_mb*1024:raise MemoryError(f'Proactive RSS guard reached: {rss} KiB > {self.rss_limit_mb} MiB')
            if s%10==0:
                mem=pathlib.Path('/sys/fs/cgroup/memory.current')
                lim=pathlib.Path('/sys/fs/cgroup/memory.max')
                if mem.exists() and lim.exists():
                    try:
                        used=int(mem.read_text());limit=int(lim.read_text())
                        if used>limit-700*1024**2:raise MemoryError(f'Proactive shared-memory guard: {used}/{limit} bytes')
                    except ValueError:pass
        self.maxa=maxa
        self.prev=prev
    def __getitem__(self,ab):
        a,b=ab
        return self.pairs[(max(a,b),min(a,b))]

def independent_reduction(a,b,c):
    """Slow original one-step Phi recurrence, for validation only, in x."""
    @functools.lru_cache(None)
    def f(a,b,c):
        if a<b or b<c:return f(*sorted((a,b,c),reverse=True))
        if b==c==0:return fmpz_poly([0]*a+[1])
        if c==0:return 2*f(a-1,b-1,1)
        return f(a+1,b-1,c-1)+f(a-1,b+1,c-1)+f(a-1,b-1,c+1)
    return f(a,b,c)

def validate_cache(cache,maxsum=16):
    count=0
    for s in range(min(maxsum,cache.maxa)+1):
        for b in range(s//2+1):
            a=s-b;p=independent_reduction(2*a,2*b,0)
            q=fmpz_poly([p[2*i] for i in range(p.degree()//2+1)])
            assert q==cache[a,b],(a,b,q,cache[a,b])
            assert all(p[i]==0 for i in range(1,p.degree()+1,2))
            count+=1
    return count

def matrix_for(d,nd,params,cache):
    small=primes(2*nd);rows=[];meta=[]
    for n,mn,g in params:
        terms={}
        for i in range(n+1):
            cni=(-1)**i*(2*n*math.comb(n+i,n-i)//(n+i))
            for j in range(i+1):
                cij=cni*math.comb(i,j)*(-4)**(i-j)
                for k,gk in enumerate(g):
                    if gk:
                        key=n-i+j+k,i
                        terms[key]=terms.get(key,0)+cij*gk
        terms={key:co for key,co in terms.items() if co}
        for m in range(mn):
            out=fmpz_poly([])
            for (a,b),co in terms.items():out+=co*cache[a+m,b]
            assert out.degree()<=n,(d,n,m,out.degree(),n)
            assert out[0]==0,(d,n,m,'constant',out[0])
            assert out[1]==12*out[2],(d,n,m,'low relation',out[1],out[2])
            assert sum(math.comb(2*i,i)*int(out[i]) for i in range(out.degree()+1))==0,(d,n,m,'central_binomial_relation')
            content=functools.reduce(math.gcd,(int(x) for x in out.coeffs()),0)
            divisor=abs(content)//smooth_strip(content,small) if content else 1
            row=[int(out[i])//divisor for i in range(3,nd+1)]
            assert all(int(out[i])%divisor==0 for i in range(out.degree()+1))
            rows.append(row);meta.append({'n':n,'m':m,'divided_small_content':str(divisor),'content_rough':str(smooth_strip(content,small))})
        print(json.dumps({'stage':'rows','d':d,'n':n,'rows_so_far':len(rows),'terms':len(terms)}),flush=True)
    return rows[::-1],meta[::-1]

def gcd_minors(rows,bound,maxminors=32,checkpoint_path=None,checkpoint_matrix_sha256=None):
    n=len(rows[0]);r=len(rows);small=primes(bound)
    if r<n:return {'status':'too_few_rows','rows':r,'columns':n}
    base=rows[:n];details=[];g=0
    if checkpoint_path is not None and checkpoint_matrix_sha256 is not None and pathlib.Path(checkpoint_path).exists():
        try:saved=json.loads(pathlib.Path(checkpoint_path).read_text())
        except (ValueError,OSError):saved={};print(json.dumps({'stage':'invalid_minor_checkpoint_ignored','path':str(checkpoint_path)}),flush=True)
        if saved.get('matrix_sha256')==checkpoint_matrix_sha256 and saved.get('bound')==bound and saved.get('rows')==r and saved.get('columns')==n:
            for j,entry in enumerate(saved['minor_certificates']):
                chosen=list(range(n))
                if j:chosen[(j-1)%n]=n+(j-1)%(r-n)
                assert entry['row_indices']==chosen and str(smooth_strip(int(entry['det']),small))==entry['rough_det']
                details.append(entry);g=math.gcd(g,int(entry['rough_det']))
            assert str(g)==saved['rough_gcd']
            print(json.dumps({'stage':'resume_exact_minors','completed_minors':len(details),'matrix_sha256':checkpoint_matrix_sha256}),flush=True)
    for j in range(len(details),maxminors):
        if g==1:break
        if j==0:chosen=list(range(n));m=base
        else:
            if r==n:break
            replace=(j-1)%n;new=n+(j-1)%(r-n)
            chosen=list(range(n));chosen[replace]=new
            m=base[:];m[replace]=rows[new]
        ts=time.time();det=int(fmpz_mat(m).det());rough=smooth_strip(det,small)
        g=math.gcd(g,rough)
        details.append({'row_indices':chosen,'det':str(det),'rough_det':str(rough),'elapsed':time.time()-ts})
        if checkpoint_path is not None:
            cp=pathlib.Path(checkpoint_path);tmp=cp.with_suffix('.tmp');tmp.write_text(json.dumps({'rows':r,'columns':n,'bound':bound,'matrix_sha256':checkpoint_matrix_sha256,'rough_gcd':str(g),'complete':g==1,'minor_certificates':details},indent=2));tmp.replace(cp)
        print(json.dumps({'stage':'minor','index':j,'det_digits':len(str(abs(det))),'rough_digits':len(str(rough)),'gcd_digits':len(str(g)),'seconds':round(time.time()-ts,4)}),flush=True)
        if g==1:break
    return {'status':'certified_no_exception_p_above_bound' if g==1 else 'remaining_rough_gcd','rough_gcd':str(g),'minor_certificates':details}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('d',type=int,nargs='+');ap.add_argument('--nd',type=int);ap.add_argument('--prefix',default='kappa0');ap.add_argument('--maxminors',type=int,default=32);ap.add_argument('--snf',action='store_true');ap.add_argument('--rss-limit-mb',type=int,default=1400)
    args=ap.parse_args();plans=[(d,*parameters(d,args.nd)) for d in args.d]
    plans.sort(key=lambda plan:max(2*n+(len(g)-1)+(mn-1) for n,mn,g in plan[2]))
    cache=ReductionCache(16,rss_limit_mb=args.rss_limit_mb);checked=validate_cache(cache)
    print(json.dumps({'stage':'cache_validated','independent_monomial_checks':checked}),flush=True)
    for d,nd,params in plans:
        maxa=max(2*n+(len(g)-1)+(mn-1) for n,mn,g in params)
        tc=time.time();cache.extend(maxa);cache_seconds=time.time()-tc
        t=time.time();rows,meta=matrix_for(d,nd,params,cache)
        obj={'d':d,'nd':nd,'kappa':0,'p_condition':f'p>{2*nd} and p not congruent +/-1 mod {2*d}','rows':rows,'row_metadata':meta,'parameters':[{'n':n,'m_count':mn,'g_T':g} for n,mn,g in params],'cache_validation_checks':checked}
        path=ROOT/f'{args.prefix}_d{d}_matrix.json';path.write_text(json.dumps(obj,separators=(',',':')))
        obj={'d':d,'nd':nd,'matrix_file':path.name,'matrix_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'shape':[len(rows),len(rows[0])],'matrix_build_seconds':time.time()-t,'additional_cache_seconds':cache_seconds,'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        obj.update(gcd_minors(rows,2*nd,args.maxminors))
        if args.snf:
            st=time.time();S=fmpz_mat(rows).snf();diag=[int(S[i,i]) for i in range(min(S.nrows(),S.ncols()))]
            obj['smith_diagonal']=[str(x) for x in diag];obj['smith_seconds']=time.time()-st
        obj['total_seconds']=time.time()-t
        (ROOT/f'{args.prefix}_d{d}_certificate.json').write_text(json.dumps(obj,indent=2))
        print(json.dumps({k:v for k,v in obj.items() if k not in ['minor_certificates','smith_diagonal']}),flush=True)

if __name__=='__main__':main()
