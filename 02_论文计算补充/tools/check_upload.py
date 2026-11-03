"""Check the exact upload file set, digests and theorem-support indexes."""
from pathlib import Path
import json,hashlib,itertools,csv
root=Path(__file__).resolve().parents[2];support=root/'02_论文计算补充'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
meta=json.loads((root/'UPLOAD_FILES.json').read_text());files=meta['files'];expected={x['path'] for x in files}
actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='UPLOAD_FILES.json' and '.git' not in p.relative_to(root).parts}
assert actual==expected,{'unexpected':sorted(actual-expected),'missing':sorted(expected-actual)}
for x in files:
 p=root/x['path'];assert p.resolve().is_relative_to(root) and p.stat().st_size==x['bytes'] and sha(p)==x['sha256'],x['path']
archives=json.loads((root/'01_五轮总材料/ARCHIVES.json').read_text());assert len(archives)==21 and len({x['round'] for x in archives})==5
assert meta['manuscript_pdf_sha256']==json.loads((support/'THEOREM_INPUTS.json').read_text())['manuscript_pdf_sha256']
rows=json.loads((support/'INPUTS.json').read_text());objects=json.loads((support/'OBJECTS.json').read_text());assert {x['sha256'] for x in rows}==set(objects)
g=json.loads((support/'indices/general_rank_59.json').read_text());assert len(g)==59 and len({x['d'] for x in g})==59 and all(x['sha256'] in objects for x in g)
b=json.loads((support/'indices/POINTWISE_BATCHES.json').read_text());assert len(b)==90 and sum(x['candidate_count'] for x in b)==1775 and all(x['manifest_sha256'] in objects for x in b)
m=json.loads((support/'indices/mixed_rules.json').read_text());pool={179,181,191,193,197,199,211,223,227,229,241};triples={tuple(sorted(x['d'])) for x in m if set(x['d'])<=pool};assert len(m)==176 and triples==set(itertools.combinations(sorted(pool),3)) and all(x['summary_sha256'] in objects for x in m)
for x in rows:assert Path(x['path']).suffix.lower() not in ['.pdf','.tex','.log','.stderr','.stdout','.a','.pyc']
end=json.loads((support/'certificates/ENDPOINT.json').read_text());assert end['P']==str(10**151) and len(end['endpoints'])==12 and len(end['gates'])==41 and all(x['pass'] for x in end['gates'])
base=json.loads((support/'certificates/FINITE_BASE_FRESH.json').read_text());assert base['status']=='PASS' and base['composite_witnesses_recomputed']==26210 and base['pointwise_source_bound_receipts']==1775
latest=json.loads((support/'certificates/FINITE_LATEST_FRESH.json').read_text());assert latest['independent_full_CRT_reconstructed'] and (latest['candidate_count'],latest['composite_witnesses'],latest['mixed_applications'],latest['direct_cases'])==(10369,9843,520,6) and latest['arithmetic']['fully_checked_dag_nodes']==18104
small=list(csv.DictReader((support/'certificates/SMALL_GRAPHS.csv').open()));assert len(small)==211 and all(x['connected']=='1' for x in small)
print(json.dumps({'status':'PASS_EXACT_UPLOAD_SET_AND_IDENTITIES','files':len(files),'round_archives':21,'unique_selected_objects':len(objects),'scope':'Artifact completeness/identity and certificate-summary consistency; not fresh regeneration of every rank determinant.'}))
