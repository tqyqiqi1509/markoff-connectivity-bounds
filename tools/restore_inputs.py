"""Fetch fixed release assets and restore a new, writable proof-input tree.

Only the Python standard library is used. Hashes are mandatory; extraction is
limited to the explicitly named content objects in PATH_MAP.json.
"""
from pathlib import Path
import argparse, hashlib, json, os, shutil, urllib.request, zipfile
from urllib.parse import quote
from bundle_schema import safe, validate

REPO = Path(__file__).resolve().parents[1]

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--asset-dir',type=Path,default=REPO/'.cache')
    ap.add_argument('--offline',action='store_true')
    a=ap.parse_args()
    if a.output.exists():raise SystemExit('Output must not already exist.')
    table=json.loads((REPO/'data/PATH_MAP.json').read_text(encoding='utf-8'))
    manifest=json.loads((REPO/'data/RELEASE_ASSETS.json').read_text(encoding='utf-8'))
    validate(table,manifest)
    a.asset_dir.mkdir(parents=True,exist_ok=True)
    for asset in manifest['assets']:
        path=safe(a.asset_dir,asset['name'])
        if not path.exists():
            if a.offline:raise SystemExit('Missing offline asset: '+str(path))
            url=f"https://github.com/{manifest['repository']}/releases/download/{quote(manifest['tag'],safe='')}/{quote(asset['name'],safe='')}"
            temp=safe(a.asset_dir,asset['name']+'.part')
            req=urllib.request.Request(url,headers={'User-Agent':'markoff-proof-supplement'})
            with urllib.request.urlopen(req,timeout=120) as response,temp.open('wb') as stream:
                shutil.copyfileobj(response,stream,1024*1024)
            if temp.stat().st_size!=asset['bytes'] or sha(temp)!=asset['sha256']:
                raise SystemExit('Downloaded asset hash mismatch: '+asset['name'])
            temp.replace(path)
        if path.stat().st_size!=asset['bytes'] or sha(path)!=asset['sha256']:
            raise SystemExit('Asset hash mismatch: '+asset['name'])
        print('verified '+asset['name'],flush=True)
    root=a.output.resolve();root.mkdir(parents=True)
    canonical={};archives={}
    try:
        for row in table['files']:
            dest=safe(root,row['path']);dest.parent.mkdir(parents=True,exist_ok=True)
            h=row['sha256']
            if h in canonical:
                # Links are internal to this writable work tree, never to the
                # repository or immutable release objects.
                try:os.link(canonical[h],dest)
                except OSError:shutil.copyfile(canonical[h],dest)
            elif row['repository_path'] is not None:
                src=safe(REPO,row['repository_path'])
                if sha(src)!=h:raise ValueError('Repository object mismatch: '+row['path'])
                shutil.copyfile(src,dest)
                canonical[h]=dest
            else:
                name=row['release_asset']
                if name not in archives:archives[name]=zipfile.ZipFile(safe(a.asset_dir,name))
                z=archives[name];member=row['member']
                if member!='objects/'+h:raise ValueError('Unexpected release member')
                with z.open(member) as source,dest.open('wb') as output:
                    shutil.copyfileobj(source,output,1024*1024)
                if sha(dest)!=h:raise ValueError('Object hash mismatch: '+row['path'])
                canonical[h]=dest
            if dest.stat().st_size!=row['bytes']:raise ValueError('Object length mismatch')
    finally:
        for z in archives.values():z.close()
    result={'status':'PASS_RESTORED_ALL_DECLARED_INPUTS','logical_paths':len(table['files']),
            'unique_objects':len(canonical),'scope':'Artifact identity and safe restoration; no mathematical certification.'}
    (root/'RESTORE_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
