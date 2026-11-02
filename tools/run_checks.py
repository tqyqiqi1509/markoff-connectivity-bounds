"""Portable entry points. Expensive stages are explicit, never automatic."""
from pathlib import Path
import argparse, os, subprocess, sys, shutil, json
from bundle_schema import safe

REPO=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=['scalar','endpoint','finite-records','resource','positive-integer','positive-nodes','mass','selected-euler','frequency','variance','relaxed-variance','small-graphs','pointwise-batch','general-modular'])
    ap.add_argument('--work-dir',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--d',type=int)
    ap.add_argument('--batch',help='Relative batch directory, from indices/POINTWISE_BATCHES.json')
    a=ap.parse_args();w=a.work_dir.resolve();out=a.output_dir.resolve()
    if not 1<=a.workers<=8:raise SystemExit('--workers must be between 1 and 8.')
    if not (w/'RESTORE_RECEIPT.json').exists():raise SystemExit('Restore the declared inputs first.')
    if out.exists():raise SystemExit('Choose a new output directory.')
    out.mkdir(parents=True);env=dict(os.environ,MARKOFF_WORKDIR=str(w),PYTHONUNBUFFERED='1')
    if sys.flags.optimize:raise SystemExit('Assertions must be enabled; do not use python -O.')
    def run(args,cwd=None):
        print('RUN',json.dumps(list(map(str,args))),flush=True)
        subprocess.run(list(map(str,args)),check=True,cwd=cwd or w,env=env)
    py=sys.executable;t=REPO/'tools';positive=w/'prior_restored/latest/positive157';pdata=positive/'complete'
    if a.stage=='scalar':run([py,t/'check_scalar_gates.py'])
    elif a.stage=='endpoint':run([py,t/'replay_continuous_budget_fraction.py','--project-root',w,'--output-dir',out/'endpoint'])
    elif a.stage=='finite-records':run([py,t/'check_finite.py','--project-root',w,'--output',out/'finite.json'])
    elif a.stage=='resource':
        for f in (t/'resource').iterdir():shutil.copyfile(f,out/f.name)
        run([py,out/'verify_resource_rational.py'],out)
    elif a.stage=='positive-integer':run([py,positive/'replay_positive157_integer.py','--directory',pdata,'--output',out/'positive.json'])
    elif a.stage=='positive-nodes':run([py,positive/'check_nodes_chebyshev.py','--directory',pdata,'--output',out/'nodes.json'])
    elif a.stage in ['mass','selected-euler','relaxed-variance']:
        name={'mass':'verify_mass_independent.py','selected-euler':'verify_selected_moment_independent.py','relaxed-variance':'verify_relaxed_budget_independent.py'}[a.stage]
        command=[py,t/name]
        if a.stage=='relaxed-variance':command+=['--workers',a.workers]
        run(command)
    elif a.stage in ['frequency','variance']:
        run([py,t/'verify_resource_gaussian_independent.py','--stage',a.stage,'--workers',a.workers])
    elif a.stage=='small-graphs':
        exe=out/('small_graphs.exe' if os.name=='nt' else 'small_graphs')
        run(['g++','-O2','-std=c++17',t/'small_graphs.cpp','-o',exe])
        with (out/'small_graphs.csv').open('w') as f:subprocess.run([str(exe)],stdout=f,check=True,env=env)
        expected=(REPO/'indices/small_graphs.csv').read_text().splitlines()
        assert (out/'small_graphs.csv').read_text().splitlines()==expected
    elif a.stage=='pointwise-batch':
        if not a.batch:raise SystemExit('--batch required')
        batches=json.loads((REPO/'indices/POINTWISE_BATCHES.json').read_text(encoding='utf-8'))
        match=[b for b in batches if b['path']==a.batch]
        if len(match)!=1:raise SystemExit('Batch must be an indexed relative path.')
        base=safe(w,a.batch);d=match[0]['d']
        run(['g++','-O2','-std=c++17',base/'independent_modular_det.cpp','-o',base/'independent_modular_det','-lgmpxx','-lgmp'])
        run([py,base/'run_manifest_group_verification.py',d,'--manifest',base/'batch_manifest.json'])
    elif a.stage=='general-modular':
        if a.d not in [139,149,151,157,163]:raise SystemExit('This entry point handles d=139,149,151,157,163; see general-rank README for earlier formats.')
        run([py,w/'agents/parallel_modular_verifier_round40/verify_universal_modular_parallel.py',a.d,'--source-dir',w/f'new_uniform/d{a.d}','--output-dir',out/f'd{a.d}','--workers',a.workers])
    fixed_results={
        'mass':['root_logs/MASS_REPLAY_PUBLIC.json'],
        'selected-euler':['root_logs/SELECTED_MOMENT_PUBLIC.json'],
        'frequency':['root_logs/RESOURCE_GAUSSIAN_FREQUENCY_ROOT_REVIEW.json'],
        'variance':['root_logs/RESOURCE_GAUSSIAN_VARIANCE_ROOT_REVIEW.json'],
        'relaxed-variance':['root_logs/RELAXED_BUDGET_VARIANCE_ROOT_REVIEW.json'],
    }
    receipts=[]
    for rel in fixed_results.get(a.stage,[]):
        source=safe(w,rel)
        if not source.is_file():raise SystemExit('Expected fresh result was not produced: '+rel)
        shutil.copyfile(source,out/source.name);receipts.append(source.name)
    (out/'RUN_STATUS.json').write_text(json.dumps({'status':'COMMAND_COMPLETED','stage':a.stage,'result_copies':receipts,'scope':'See README for exactly what this stage checks. Analytic stages and pointwise batches can write fixed result names in the restored writable work tree; immutable repository/release inputs are not modified.'},indent=2)+'\n')

if __name__=='__main__':main()
