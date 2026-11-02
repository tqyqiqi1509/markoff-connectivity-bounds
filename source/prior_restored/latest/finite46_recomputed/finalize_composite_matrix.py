#!/usr/bin/env python3
"""Assemble resumable phi_stream rows and certify maximal minors."""
import argparse, functools, hashlib, json, math, pathlib, resource, signal, subprocess, sys, time
from martin_kappa0 import primes,smooth_strip,gcd_minors,factor
from composite_parameters import parameters
sys.set_int_max_str_digits(0)
ROOT=pathlib.Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('d',type=int);ap.add_argument('raw');ap.add_argument('--prefix',default='composite');ap.add_argument('--reference');ap.add_argument('--maxminors',type=int,default=32);ap.add_argument('--method',choices=['auto','minors','codimension_one'],default='auto');args=ap.parse_args()
ts=time.time();rawpath=pathlib.Path(args.raw);allraw=[json.loads(line) for line in rawpath.read_text().splitlines() if line.strip()];assert allraw
nd=allraw[0]['nd'];nd,params=parameters(args.d,nd)
cpp_plan=[json.loads(line) for line in pathlib.Path(str(rawpath)+'.parameters.jsonl').read_text().splitlines()]
assert [(p['n'],p['m_count'],list(map(int,p['g_T']))) for p in cpp_plan]==params
expected={(n,m) for n,mn,g in params for m in range(mn)}
indexed={}
for r in allraw:
    assert r['d']==args.d and r['nd']==nd
    key=r['n'],r['m'];assert key in expected and key not in indexed;indexed[key]=r
assert set(indexed)==expected,{'missing':sorted(expected-set(indexed))}
rows=[];meta=[];ps=primes(2*nd)
for key in sorted(expected,reverse=True):
    r=indexed[key];co=list(map(int,r['coefficients']));assert len(co)==r['n']+1
    assert co[0]==0 and co[1]==12*co[2]
    assert sum(math.comb(2*i,i)*c for i,c in enumerate(co))==0
    content=functools.reduce(math.gcd,co,0);rough=smooth_strip(content,ps);div=abs(content)//rough if content else 1
    row=[(co[i] if i<len(co) else 0)//div for i in range(3,nd+1)]
    assert all(c%div==0 for c in co)
    rows.append(row);meta.append({'n':r['n'],'m':r['m'],'divided_small_content':str(div),'content_rough':str(rough)})
comparison=None
if args.reference:
    ref=json.loads(pathlib.Path(args.reference).read_text());assert rows==ref['rows'];assert meta==ref['row_metadata'];comparison={'file':args.reference,'all_rows_and_content_normalizations_equal':True}
obj={'d':args.d,'nd':nd,'kappa':0,'general_d_theorem_status':'conditional on parent general-d bootstrap audit','filter_audit_file':f'composite_d{args.d}_parameters_audit.json','p_condition':f'p>{2*nd} and p not congruent +/-1 mod {2*args.d}','rows':rows,'row_metadata':meta,'parameters':[{'n':n,'m_count':mn,'g_T':g} for n,mn,g in params],
     'engine':'C++ GMP scalar coefficient propagation through descending total degree layers','engine_source_sha256':hashlib.sha256((ROOT/'phi_stream_composite.cpp').read_bytes()).hexdigest(),
     'raw_rows_file':rawpath.name,'raw_rows_sha256':hashlib.sha256(rawpath.read_bytes()).hexdigest(),'reference_comparison':comparison}
mp=ROOT/f'{args.prefix}_d{args.d}_matrix.json';mp.write_text(json.dumps(obj,separators=(',',':')))
cert={'d':args.d,'nd':nd,'matrix_file':mp.name,'matrix_sha256':hashlib.sha256(mp.read_bytes()).hexdigest(),'shape':[len(rows),len(rows[0])],
      'matrix_build_seconds':time.time()-ts,'engine':obj['engine'],'reference_comparison':comparison,
      'row_compute_seconds_sum':sum(r['seconds'] for r in allraw),'nonzero_states_total':sum(r['nonzero_states'] for r in allraw),'integer_additions_total':sum(r['additions'] for r in allraw)}
method=args.method
if method=='auto':method='codimension_one' if args.d>=67 and factor(args.d)==[(args.d,1)] and args.prefix=='stream' else 'minors'
if method=='codimension_one':
    assert args.prefix=='stream' and factor(args.d)==[(args.d,1)] and args.d>=5
    child=subprocess.Popen([sys.executable,str(ROOT/'codimension_one_certificate.py'),str(args.d)],cwd=ROOT)
    try:
        rc=child.wait()
    except KeyboardInterrupt:
        if child.poll() is None:
            child.send_signal(signal.SIGINT)
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:child.kill();child.wait()
        raise
    if rc:raise subprocess.CalledProcessError(rc,child.args)
    cp=ROOT/f'codim_d{args.d}_certificate.json';cc=json.loads(cp.read_text())
    assert cc['source_sha256']==cert['matrix_sha256']
    if cc['status']=='certified_no_exception_p_above_bound':
        cert.update({'status':cc['status'],'rough_gcd':cc['rough_gcd'],'certificate_method':'codimension_one','codimension_one_certificate_file':cp.name,'codimension_one_certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest()})
    elif cc['status']=='lower_block_obstruction' and args.method=='auto':
        print(json.dumps({'stage':'codimension_one_inconclusive','d':args.d,'fallback':'full maximal minors'}),flush=True);method='minors'
    else:cert.update({'status':cc['status'],'rough_gcd':cc.get('rough_gcd'),'certificate_method':'codimension_one','codimension_one_certificate_file':cp.name,'codimension_one_certificate_sha256':hashlib.sha256(cp.read_bytes()).hexdigest()})
if method=='minors':cert.update(gcd_minors(rows,2*nd,args.maxminors,ROOT/f'{args.prefix}_d{args.d}_minors_checkpoint.json',cert['matrix_sha256']));cert['certificate_method']='maximal_minors'
cert['finalization_seconds']=time.time()-ts;cert['maxrss_kib']=max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
cp=ROOT/f'{args.prefix}_d{args.d}_certificate.json';cp.write_text(json.dumps(cert,indent=2))
print(json.dumps({k:v for k,v in cert.items() if k!='minor_certificates'}),flush=True)
