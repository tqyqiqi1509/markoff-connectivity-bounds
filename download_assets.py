"""Download the two-layer release attachments and verify exact SHA-256 values."""
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parent
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--layer',choices=['1','2','all'],default='all');ap.add_argument('--list',action='store_true');args=ap.parse_args()
 assets=json.loads((ROOT/'RELEASE_ASSETS.json').read_text())['assets']
 for a in assets:
  rel=a['path'];assert rel.startswith(('01_','02_')) and rel.endswith('.zip');p=(ROOT/rel).resolve();assert p.is_relative_to(ROOT)
  if args.layer!='all' and not rel.startswith('0'+args.layer+'_'):continue
  if args.list:print(rel,a['bytes']);continue
  if p.exists():
   if p.stat().st_size!=a['bytes'] or digest(p)!=a['sha256']:raise RuntimeError('Existing file differs: '+rel)
   print('VERIFIED',rel,flush=True);continue
  p.parent.mkdir(parents=True,exist_ok=True);partial=p.with_suffix(p.suffix+'.partial')
  for attempt in range(5):
   try:
    start=partial.stat().st_size if partial.exists() else 0
    if start==a['bytes'] and digest(partial)==a['sha256']:partial.rename(p);break
    if start>=a['bytes']:raise RuntimeError('Partial file differs; inspect '+str(partial))
    headers={'User-Agent':'Markoff-proof-materials-v3'}
    if start:headers['Range']='bytes='+str(start)+'-'
    with urlopen(Request(a['url'],headers=headers),timeout=120) as response:
     resumed=response.status==206 and start>0;count=start if resumed else 0
     with partial.open('ab' if resumed else 'wb') as f:
      while True:
       block=response.read(4*1024*1024)
       if not block:break
       count+=len(block)
       if count>a['bytes']:raise RuntimeError('Unexpected download size: '+rel)
       f.write(block)
    if count!=a['bytes'] or digest(partial)!=a['sha256']:raise RuntimeError('Incomplete or mismatched download: '+rel)
    partial.rename(p);print('VERIFIED',rel,flush=True);break
   except (HTTPError,URLError,TimeoutError,OSError) as error:
    if attempt==4:raise
    print('RETRY',rel,type(error).__name__,flush=True);time.sleep(3*(attempt+1))
if __name__=='__main__':main()
