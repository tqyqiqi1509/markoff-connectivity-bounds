from pathlib import Path
import sys,json,hashlib,math,time
from flint import fmpz_mat
sys.set_int_max_str_digits(0)
R=Path(__file__).resolve().parent.parent;D=R/'objects';out=R/'certificates';out.mkdir(exist_ok=True)
def read(h):
 p=D/h;b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==h;return json.loads(b)
s=read('df21615262643c327ef68a248343276094c8a9b5e0b7cb347bb14c07e942cf8f')
rows={};source_receipts=[]
for r in s['source_bindings']:
 b=(D/r['sha256']).read_bytes();assert hashlib.sha256(b).hexdigest()==r['sha256']
 for line in b.splitlines():
  o=json.loads(line);key=(o['d'],o['n'],o['m']);assert key not in rows;rows[key]=list(map(int,o['coefficients']))
 source_receipts.append(r)
assert len(rows)==s['source_rows']==194
A=[]
for key,norm in zip(s['row_keys'],s['row_normalizers']):
 v=rows[tuple(key)];assert v[0]==0 and v[1]==12*v[2];v=v+[0]*(s['N']+1-len(v));assert len(v)==s['N']+1
 assert sum(value*math.comb(2*j,j) for j,value in enumerate(v))==0
 g=math.gcd(*v);norm=int(norm);assert g==norm and all(x%norm==0 for x in v)
 rough=norm
 for q in range(2,s['strict_p_gate']+1):
  while rough%q==0:rough//=q
 assert rough==1
 A.append([x//norm for x in v[3:]])
assert len(A)==194 and all(len(r)==184 for r in A)
results=[];dets=[]
for i,r in enumerate(s['minor_certificates']):
 cert=read(r['sha256']);ix=cert['selected_row_indices'];assert len(ix)==len(set(ix))==184
 M=[A[j] for j in ix];H=1
 for row in M:
  n=sum(x*x for x in row);b=math.isqrt(n);H*=b+(b*b<n)
 start=time.monotonic();det=int(fmpz_mat(M).det());assert det==int(cert['integer_determinant']) and abs(det)<=H
 dets.append(det);results.append({'minor':i,'dimension':184,'determinant_sha256':hashlib.sha256(str(det).encode()).hexdigest(),'determinant_bits':abs(det).bit_length(),'seconds':time.monotonic()-start,'direct_integer_determinant_matches':True})
 print('EXACT_MINOR_PASS',i,flush=True)
g=math.gcd(*dets);assert g==int(s['gcd']);factors=[];r=g
for q in range(2,s['strict_p_gate']+1):
 a=0
 while r%q==0:r//=q;a+=1
 if a:factors.append([q,a])
assert r==1
receipt={'status':'PASS_NEW_EXACT_SOURCE_MATRIX_TO_TWO_INTEGER_MINORS_CONNECTION','case':101,'N':186,'gate':373,'anchor':62,'required_legal_d':[62,179],'source_rows':194,'columns':184,'source_bindings':source_receipts,'minors':results,'gcd_support':factors,'rough_gcd_above_gate':1,'scope':'Fresh complete row-to-matrix connection and direct fmpz determinants. Legal raw-source construction follows the supplied source recurrence/authentication; those raw source coefficients were not regenerated here. This new certificate replaces a missing historical connection note, whose expected bytes remain unrecovered.'}
(out/'CASE101_RECONNECTION.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
print(receipt['status'],flush=True)
