#!/usr/bin/env python3
"""General-d Chebyshev filters, independently checked in cyclotomic fields.

Source: Martin section8.1 defines the filter for arbitrary d; the root
worker is separately auditing the general-d bootstrap theorem.  This file
checks source-filter arithmetic only, not connectivity.
"""
import argparse,datetime,hashlib,json,math,pathlib,sys,time
sys.set_int_max_str_digits(0)
from martin_kappa0 import ROOT,cheb_even,factor,val,parameters as prime_power_parameters
from flint import fmpz_poly

def parameters(d,nd):
    if d<4 or nd<3*d or nd%d:raise ValueError('require d>=4 and nd>=3d a multiple of d')
    result=[]
    for n in range(d,nd+1,d):
        G=fmpz_poly([1])
        for ell,a in factor(d):
            r=n//ell**(val(n,ell)-a+1)
            U=fmpz_poly(cheb_even(r-1))
            G=G*U//G.gcd(U)
        H=fmpz_poly(cheb_even(n-1));quotient,remainder=divmod(H,G)
        assert not remainder and G[G.degree()]==1
        count=n//2-G.degree();assert count==quotient.degree()>0
        result.append((n,count,[int(c) for c in G.coeffs()]))
    return nd,result

def audit(d,nd):
    st=time.time();_,params=parameters(d,nd);records=[]
    for n,count,co in params:
        z=fmpz_poly([0,1]);modulus=fmpz_poly.cyclotomic(2*n);G=fmpz_poly(co);H=fmpz_poly(cheb_even(n-1))
        assert H.gcd(H.derivative()).degree()==0
        powers=[fmpz_poly([1])]
        for i in range(1,2*n):powers.append(powers[-1]*z%modulus)
        selected=[];excluded=[];seen=set()
        for i in range(1,n//2+1):
            beta=(powers[(2*i)%(2*n)]+powers[(-2*i)%(2*n)]+2)%modulus
            key=tuple(map(int,beta.coeffs()));assert key not in seen;seen.add(key)
            def evaluate(P):
                out=fmpz_poly([])
                for c in reversed(P.coeffs()):out=(out*beta+c)%modulus
                return out
            assert not evaluate(H)
            ordinary_order=2*n//math.gcd(i,2*n)
            even_order=math.lcm(ordinary_order,2)
            allowed=even_order%(2*d)==0
            assert (not evaluate(G))== (not allowed),(d,n,i)
            (selected if allowed else excluded).append({'representative_exponent':i,'ordinary_order':ordinary_order,'even_order':even_order})
        assert len(selected)==count and len(excluded)==G.degree()
        records.append({'n':n,'m_count':count,'g_T':co,'cyclotomic_modulus_degree':modulus.degree(),'selected_beta_classes':selected,'excluded_beta_classes':excluded,'exact_filter_zero_set_verified':True})
    result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'d':d,'nd':nd,'matrix_shape':[sum(m for n,m,g in params),nd-2],'source':'Martin section8.1 arbitrary-d filter; no connectivity theorem is asserted by this audit.','filter_algorithm':'Python monic polynomial lcm, independently compared with exact evaluation at all squared traces in Q[zeta_(2n)]','records':records,'seconds':time.time()-st}
    (ROOT/f'composite_d{d}_parameters_audit.json').write_text(json.dumps(result,indent=2));return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('plans',nargs='+',help='d:nd');args=ap.parse_args()
    old=[]
    for d in [4,5,7,8,9,11,13,16,17,19]:
        nd,p=prime_power_parameters(d);_,q=parameters(d,nd);assert p==q;old.append(d)
    (ROOT/'composite_prime_power_regression.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'old_prime_power_parameters_exactly_equal':old},indent=2))
    for plan in args.plans:
        d,nd=map(int,plan.split(':'));r=audit(d,nd)
        print(json.dumps({k:v for k,v in r.items() if k!='records'}),flush=True)

if __name__=='__main__':main()
