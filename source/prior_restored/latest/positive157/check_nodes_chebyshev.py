#!/usr/bin/env python3
"""Independent full pointwise verification by real Chebyshev/geometric formula.
Producer uses direct acb polynomial addition with root powers. This checker
uses exp(log) constants and |sum z^j|^2=(1+r^(2e+2)-2r^(e+1)T_(e+1)(cos t))/(1+r^2-2r cos t).
Both computations use the Arb library; this is not an independent library claim.
Every retained endpoint/tail is checked against its saved integer upper bound.
"""
import argparse,gzip,hashlib,json,struct,time,sys
from pathlib import Path
from flint import arb,ctx
import flint
sys.set_int_max_str_digits(0)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def constants(rows,prec):
 ctx.prec=prec
 out=[]
 for x in rows[1:]:
  q=x['q'];log=arb(q).log();r=(-log/17).exp();s=(-log/18).exp()
  B=arb(1);sp=arb(1);co=[];powers=[]
  for e in range(1,x['first_discarded_exponent']):
   sp*=s;B+=(arb(q-1)/q)*sp
   co.append((-arb(829)*e*log/8100+arb(2)*B.log()/15).exp())
   powers.append((r**(e+1),r**(2*e+2)))
  out.append({'q':q,'log':log,'r':r,'r2':r*r,'coeff':co,'powers':powers})
 return out
def node_upper(data,i,N,density,rank,prec):
 ctx.prec=prec
 t=arb.pi()*i/N if rank else arb(i)/density
 local=[]
 for row in data:
  c=(t if rank else t*row['log']).cos()
  den=1+row['r2']-2*row['r']*c;assert den.lower()>0
  cm1=arb(1);cm=c;g=arb(1)
  for coeff,(rn,rn2) in zip(row['coeff'],row['powers']):
   cn=2*c*cm-cm1;cm1,cm=cm,cn
   num=1+rn2-2*rn*cn;assert num.lower()>0
   val=coeff*(num/den).sqrt();g=max(g,val.upper())
  local.append(g)
 result=[arb(1)]*(len(local)+1);prod=arb(1)
 for k in range(len(local)-1,-1,-1):
  prod=(prod*local[k]).upper();result[k]=prod
 return result
def run(folder,mode,p,metadata):
 info=json.loads((folder/(mode+'.json')).read_text());path=folder/info['node_file']
 assert sha(path)==info['node_gzip_sha256']
 N=info['panels'];columns=info['columns'];NS=1<<info['node_scale_bits'];codec=struct.Struct('<'+'Q'*columns)
 data=constants(p['active_primes'],384);prec=256;fallback=[];rawhash=hashlib.sha256();start=time.monotonic()
 with gzip.open(path,'rb') as stream:
  for i in range(N+1):
   raw=stream.read(codec.size);assert len(raw)==codec.size;rawhash.update(raw);values=codec.unpack(raw)
   bounds=node_upper(data,i,N,metadata['density'],mode=='rank',prec)
   assert len(bounds)==columns==len(p['active_primes'])
   failed=[j for j,b in enumerate(bounds) if not b<=arb(values[j])/NS]
   if failed:
    more=constants(p['active_primes'],768);bounds=node_upper(more,i,N,metadata['density'],mode=='rank',512)
    assert all(b<=arb(values[j])/NS for j,b in enumerate(bounds)),(mode,i,failed)
    fallback.append({'node':i,'columns':failed,'second_precision':512})
   if (i+1)%16384==0 or i==N:
    print(json.dumps({'mode':mode,'nodes_checked':i+1,'total_nodes':N+1,'seconds':time.monotonic()-start}),flush=True)
  assert stream.read(1)==b''
 assert rawhash.hexdigest()==info['node_raw_sha256']
 return {'mode':mode,'every_pointwise_tail_upper_verified':True,'nodes':N+1,'tail_values':(N+1)*columns,
  'normal_precision':prec,'constant_precision':384,'high_precision_fallbacks':fallback,
  'node_raw_sha256':rawhash.hexdigest(),'node_gzip_sha256':sha(path),'seconds_display_only':time.monotonic()-start}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--directory',required=True);ap.add_argument('--output',required=True)
 args=ap.parse_args();folder=Path(args.directory);header=json.loads((folder/'parameters.json').read_text())
 checks=[run(folder,mode,header['parameters'],header['setup']) for mode in ['rank','mellin']]
 result={'status':'PASS','all_saved_transcendental_nodes_recomputed':True,
  'independent_formula':'real Chebyshev closed geometric modulus squared, exp(log) local coefficients',
  'shared_numeric_library':'python-flint Arb, not independent libraries','flint_version':flint.__version__,
  'checks':checks,'parameters_sha256':sha(folder/'parameters.json'),'script_sha256':sha(Path(__file__))}
 Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
