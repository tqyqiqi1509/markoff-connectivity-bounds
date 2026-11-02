"""Verify every mapped repository/release input, without running producers."""
from pathlib import Path
import argparse, hashlib, json, zipfile, ast, sys
from bundle_schema import safe, validate

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    if sys.flags.optimize:raise SystemExit('Assertions must be enabled; do not use python -O.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--asset-dir',type=Path)
    ap.add_argument('--work-dir',type=Path)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();repo=Path(__file__).resolve().parents[1]
    m=json.loads((repo/'data/PATH_MAP.json').read_text(encoding='utf-8'))
    assets=json.loads((repo/'data/RELEASE_ASSETS.json').read_text(encoding='utf-8'))
    validate(m,assets)
    checked={};zs={};binary_count=0
    try:
        if a.asset_dir:
            for x in assets['assets']:
                p=safe(a.asset_dir,x['name'])
                assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256'],x['name']
                zs[x['name']]=zipfile.ZipFile(p)
                assert zs[x['name']].testzip() is None
        for r in m['files']:
            key=(r['repository_path'],r['sha256'])
            if key not in checked:
                if r['repository_path']:
                    p=safe(repo,r['repository_path']);assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['path']
                    if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
                elif a.asset_dir:
                    z=zs[r['release_asset']]
                    assert z.getinfo(r['member']).file_size==r['bytes']
                    with z.open(r['member']) as f:assert hashlib.file_digest(f,'sha256').hexdigest()==r['sha256']
                elif not a.work_dir:continue
                if r['repository_path'] or a.asset_dir:checked[key]=True
            if a.work_dir:
                p=safe(a.work_dir,r['path']);assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['path']
        general=json.loads((repo/'indices/general_rank_inputs.json').read_text(encoding='utf-8'))
        bypath={r['path']:r for r in m['files']}
        assert len(general)==57
        for r in general:assert bypath[r['certificate_path']]['sha256']==r['sha256']
    finally:
        for z in zs.values():z.close()
    result={'status':'PASS_ARTIFACT_IDENTITIES' if a.asset_dir or a.work_dir else 'PASS_REPOSITORY_IDENTITIES_ONLY','logical_paths':len(m['files']),
            'checked_storage_objects':len(checked),'all_release_assets_checked':bool(a.asset_dir),
            'all_restored_input_paths_checked':bool(a.work_dir),
            'all_57_general_certificate_paths_resolved':True,'scope':'Identity, syntax and declared coverage only; no fresh large determinant computations.'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
