"""Expand the proof-data volumes with path, size and SHA-256 checks."""
from pathlib import Path
import json,zipfile,hashlib
r=Path(__file__).resolve().parent.parent;meta=json.loads((r/'DATA_VOLUMES.json').read_text());objects=json.loads((r/'OBJECTS.json').read_text());seen=set()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
for v in meta['volumes']:
 p=r/v['file'];assert sha(p)==v['sha256']
 with zipfile.ZipFile(p) as z:
  for i in z.infolist():
   name=i.filename;assert name.startswith('objects/') and name.count('/')==1
   h=name.split('/')[1];assert h in objects and h not in seen and i.file_size==objects[h]['bytes'];seen.add(h)
   dst=r/name;assert dst.resolve().is_relative_to(r);dst.parent.mkdir(exist_ok=True)
   if dst.exists():assert sha(dst)==h;continue
   dig=hashlib.sha256()
   with z.open(i) as f,dst.open('xb') as g:
    for b in iter(lambda:f.read(8*1024*1024),b''):dig.update(b);g.write(b)
   assert dig.hexdigest()==h
assert seen==set(objects)
print('PASS_DATA_EXPANSION',len(seen))
