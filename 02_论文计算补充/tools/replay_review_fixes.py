"""Replay the review fixes: dependency binding and the three middle arithmetic intervals.

This does not regenerate the 59 general or 176 mixed determinant certificates.
Run after prepare_data.py and materialize.py. Optional matrix builds reconstruct
only d=131 and d=137; optional case9002 needs Linux and python-flint==0.9.0.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--build-131-137',action='store_true')
    ap.add_argument('--case9002',action='store_true')
    a=ap.parse_args();root=a.project_root.resolve();out=a.output.resolve()
    assert root.is_dir() and not out.exists()
    index_path=Path(__file__).resolve().parent.parent/'indices/REVIEW_DEPENDENCIES.json'
    index=json.loads(index_path.read_text(encoding='utf-8'))
    for x in index['supplied']:
        p=root/x['path'];assert p.is_file() and p.stat().st_size==x['bytes']
        assert digest(p)==x['sha256'],x['path']
    for x in index['parameter_sources']:
        source=root/x['source'];parameter=Path(str(source)+'.parameters.json')
        assert digest(source)==x['source_sha256'] and digest(parameter)==x['parameters_sha256']
        meta=json.loads(parameter.read_text());rows=[json.loads(line) for line in source.open()]
        assert (meta['d'],meta['n'],meta['m_count'])==(x['d'],x['n'],x['rows'])
        assert len(rows)==x['rows']
        assert all((r['d'],r['n'],r['m'])==(x['d'],x['n'],i) for i,r in enumerate(rows))
    out.mkdir(parents=True)
    spec=index['middle_intervals'];summary=[]
    wrapper=root/'agents/finite_modular/replay_named_arithmetic.py'
    assert digest(wrapper)==index['middle_wrapper_sha256']
    engine=root/'agents/window_independent_review/independent_window_replay.py'
    assert digest(engine)==index['middle_engine_sha256']
    for x in spec:
        target=out/x['id']
        cmd=[sys.executable,str(wrapper),'--source',str(root/x['source']),
             '--window',str(root/x['window']),'--dag',str(root/x['dag']),'--out',str(target)]
        for field in ('source','window','dag'):
            assert digest(root/x[field])==x[field+'_sha256']
        subprocess.run(cmd,check=True)
        receipt=json.loads((target/'REVIEW.json').read_text())
        assert receipt['direct_passed']==x['direct'] and receipt['window_passed']==x['window_count']
        summary.append({'id':x['id'],'direct':receipt['direct_passed'],'window':receipt['window_passed'],
                        'dag_nodes':receipt['independently_checked_dag_nodes'],
                        'receipt_sha256':digest(target/'REVIEW.json')})
    assert sum(x['direct'] for x in summary)==397 and sum(x['window'] for x in summary)==57
    builds=[]
    if a.build_131_137:
        for x in index['generated']:
            path=root/x['path'];assert not path.exists(),('generated output exists',str(path))
            script=root/x['builder'];subprocess.run([sys.executable,str(script),str(x['d'])],check=True)
            assert digest(path)==x['sha256']
            builds.append({'d':x['d'],'sha256':digest(path)})
    mixed=None
    if a.case9002:
        # The frozen checker binds output paths relative to its project root.
        assert out.is_relative_to(root),'For --case9002 choose --output inside --project-root.'
        cmd=[sys.executable,str(root/'agents/schur_integer_audit/independent_compressed_mixed_check.py'),
             '--plan',str(root/'computations/mixed_schur/PLAN_STAR169.json'),'--case','9002',
             '--producer',str(root/'computations/mixed_schur/case9002'),'--out',str(out/'case9002')]
        subprocess.run(cmd,check=True)
        mixed=json.loads((out/'case9002/COMPRESSED_MIXED_REVIEW.json').read_text())['status']
    result={'status':'PASS_REVIEW_DEPENDENCIES_AND_MIDDLE_ARITHMETIC','supplied_paths_checked':len(index['supplied']),
            'middle_intervals':summary,'matrix_builds':builds,'case9002':mixed,
            'scope':'Exact supplied-file identities and 397 direct/57 window arithmetic replays with prime-DAG checks. Optional indicated builds and one mixed case only. No all-rank regeneration or formal proof certification is claimed.'}
    (out/'REVIEW_FIXES.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
