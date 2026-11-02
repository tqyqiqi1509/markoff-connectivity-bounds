#!/usr/bin/env python3
"""Use an already certified base determinant and one exact rational solve.

For the base square A and u=A^{-1}e_i, replacing row i by b gives
det(A_i(b))=det(A)*(b.u). A single cofactor column can certify many replacement
rows without recomputing a second full determinant. The output is separate
from production certificates until the caller reviews its exact identities.
"""
import argparse, datetime, functools, hashlib, json, math, pathlib, sys, time
sys.set_int_max_str_digits(0)
ROOT=pathlib.Path(__file__).resolve().parents[1]
COMPUTE=ROOT/'agents/computational_mechanisms'
sys.path.insert(0,str(COMPUTE))
from martin_kappa0 import primes, smooth_strip
from flint import fmpz_mat, ctx

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def atomic(path,obj):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(obj,separators=(',',':'))+'\n')
    temporary.replace(path)
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('source',type=pathlib.Path)
    ap.add_argument('--base-minors',type=pathlib.Path,required=True)
    ap.add_argument('--output',type=pathlib.Path,required=True)
    ap.add_argument('--solve-checkpoint',type=pathlib.Path,required=True)
    ap.add_argument('--column',type=int,default=0)
    ap.add_argument('--threads',type=int,default=2)
    args=ap.parse_args(); start=time.monotonic(); ctx.threads=args.threads
    source_sha=sha(args.source); src=json.loads(args.source.read_text())
    rows=src['rows']; n=len(rows[0]); i=args.column
    assert 0<=i<n and len(rows)>n and all(len(r)==n for r in rows)
    base_sha=sha(args.base_minors); base=json.loads(args.base_minors.read_text())
    assert base['matrix_sha256']==source_sha
    base_entry=base['minor_certificates'][0]
    assert base_entry['row_indices']==list(range(n))
    determinant=int(base_entry['det']); assert determinant
    binding={'method_version':2,'equation':'A Q = D e_i',
             'source_sha256':source_sha,'base_minors_sha256_at_start':base_sha,
             'base_determinant_sha256':hashlib.sha256(str(determinant).encode()).hexdigest(),
             'source_file':str(args.source),'n':n,'column':i,'d':src['d']}
    A=fmpz_mat(rows[:n])
    if args.solve_checkpoint.exists():
        saved=json.loads(args.solve_checkpoint.read_text())
        assert saved['binding']==binding
        D=int(saved['denominator']); Q=list(map(int,saved['numerators']))
        print(json.dumps({'stage':'reuse_saved_cofactor_solve','utc':stamp()}),flush=True)
    else:
        rhs=fmpz_mat([[int(j==i)] for j in range(n)])
        print(json.dumps({'stage':'one_rhs_rational_solve','utc':stamp(),'n':n,'threads':args.threads}),flush=True)
        solved=A.solve(rhs)
        D=1
        for j in range(n): D=math.lcm(D,int(solved[j,0].denominator))
        Q=[int(solved[j,0].numerator)*(D//int(solved[j,0].denominator)) for j in range(n)]
        assert D>0 and functools.reduce(math.gcd,Q,D)==1
        saved={'status':'candidate_equation_not_yet_verified','binding':binding,
               'denominator':str(D),'numerators':list(map(str,Q)),
               'solve_elapsed_seconds':time.monotonic()-start,'created_utc':stamp()}
        atomic(args.solve_checkpoint,saved)
        print(json.dumps({'stage':'candidate_QD_saved','utc':stamp(),'denominator_digits':len(str(D))}),flush=True)
    assert D>0 and len(Q)==n and functools.reduce(math.gcd,Q,D)==1
    assert A*fmpz_mat([[q] for q in Q])==fmpz_mat([[D*int(j==i)] for j in range(n)])
    assert determinant%D==0
    saved['status']='exact_cofactor_equation_verified'
    saved['verified_utc']=stamp(); atomic(args.solve_checkpoint,saved)
    small=primes(2*src['nd']); gcd_value=smooth_strip(determinant,small)
    replacements=[]
    for j in range(n,len(rows)):
        dot=sum(a*b for a,b in zip(rows[j],Q))
        value=(determinant//D)*dot
        rough=smooth_strip(value,small)
        gcd_value=math.gcd(gcd_value,rough)
        chosen=list(range(n)); chosen[i]=j
        replacements.append({'row_indices':chosen,'det':str(value),'rough_det':str(rough),
                             'derivation':'base determinant times row dot exact cofactor column'})
        if gcd_value==1: break
    assert sha(args.source)==source_sha
    # A concurrent process must not alter the imported base certificate silently.
    assert sha(args.base_minors)==base_sha
    result={'status':'certified_no_exception_p_above_bound' if gcd_value==1 else 'inconclusive_rough_gcd',
            'binding':binding,'bound':2*src['nd'],'nd':src['nd'],
            'base_determinant_imported_not_recomputed':True,
            'base_determinant':str(determinant),'denominator':str(D),
            'cofactor_checkpoint_file':str(args.solve_checkpoint),
            'cofactor_checkpoint_sha256':sha(args.solve_checkpoint),
            'exact_equation':'A Q = D e_i','rough_gcd':str(gcd_value),
            'replacement_minor_certificates':replacements,'elapsed_seconds':time.monotonic()-start,
            'finished_utc':stamp()}
    atomic(args.output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['base_determinant','denominator','replacement_minor_certificates']}),flush=True)

if __name__=='__main__': main()
