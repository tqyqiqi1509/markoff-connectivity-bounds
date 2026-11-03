"""Restore the selected logical paths; conflicting historical versions stay explicit."""
from pathlib import Path,PurePosixPath
import json,sys,os,shutil,hashlib
r=Path(__file__).resolve().parent.parent;dst=Path(sys.argv[1]).resolve();dst.mkdir(parents=True,exist_ok=True);conflicts=[];made=0
for x in json.loads((r/'INPUTS.json').read_text()):
 n=x['path'];parts=PurePosixPath(n).parts
 if n.startswith('/') or ':' in n or '..' in parts:raise ValueError(n)
 p=dst.joinpath(*parts);assert p.resolve().is_relative_to(dst);src=r/'objects'/x['sha256']
 if p.exists():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=x['sha256']:conflicts.append(x)
 else:
  p.parent.mkdir(parents=True,exist_ok=True)
  try:os.link(src,p)
  except OSError:shutil.copyfile(src,p)
  made+=1
(dst/'VERSION_CONFLICTS.json').write_text(json.dumps(conflicts,indent=2))
print(json.dumps({'restored_paths':made,'historical_conflicting_versions':len(conflicts),'selection':'first-listed; every alternate remains in objects by SHA-256'}))
