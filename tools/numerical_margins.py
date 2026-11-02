#!/usr/bin/env python3
"""Supplementary margin calculation; Python standard library only.

Decimal approximations propose rational brackets; exact Fraction power
comparisons decide all twelve endpoint brackets. The original integer gates
are checked exactly. Scalar-margin decimals are explanatory displays only.
This tool reads the endpoint receipt and never changes manuscript files.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext, ROUND_FLOOR
import json,hashlib
import argparse,sys
if not __debug__:
 raise SystemExit('Assertions must remain enabled; do not run Python with -O.')
ap=argparse.ArgumentParser(description='Exact rational margins for the twelve endpoint rows and exact scalar gates of the frozen Markoff paper.')
ap.add_argument('--endpoint-receipt',required=True,type=Path,help='FRACTION_CONTINUOUS_BUDGET_REVIEW.json from the supplement')
ap.add_argument('--output-dir',required=True,type=Path,help='A new directory; existing output directories are refused')
args=ap.parse_args()
source=args.endpoint_receipt.resolve()
O=args.output_dir.resolve()
if O.exists():raise SystemExit('Choose a new output directory; no overwrite.')
if not source.is_file():raise SystemExit('Endpoint receipt is not a regular input file.')
O.mkdir(parents=True)
j=json.loads(source.read_text(encoding='utf-8'))
def q(x):return F(int(x['numerator']),int(x['denominator']))
def dd(x):return D(x.numerator)/D(x.denominator)
assert j['status']=='PASS_FIXED_CONTINUOUS_BUDGET_ENDPOINT_PURE_INTEGER_FRACTION'
assert int(j['P'])==8*10**152 and int(j['C'])==28885
assert len(j['endpoint_rows'])==12
rows=[]
with localcontext() as c:
 c.prec=75
 for row in j['endpoint_rows']:
  u,w,t=[q(row[k]) for k in ['U_upper','W_upper','trace_upper']]
  val=((D(1510)*dd(u).ln()+D(377)*dd(w).ln()+D(252)*dd(t).ln())/D(1887)).exp()
  k=int((val*10**6).to_integral_value(rounding=ROUND_FLOOR))
  lo,hi=F(k,10**6),F(k+1,10**6)
  power=u**1510*w**377*t**252
  assert lo**1887<=power<hi**1887
  lower,upper=F(28885)-hi,F(28885)-lo
  assert lower>0
  rows.append({'region':row['region'],'point':row['point'],'F_upper_enclosure':[str(lo),str(hi)],'margin_from_28885_enclosure':[str(lower),str(upper)],'margin_lower_decimal':str(dd(lower)),'relative_margin_lower':str(dd(lower)/D(28885)),'exact_rational_power_bracket_checked':True})
 P=8*10**152
 ordinary=(D(21718).ln()+D(2541935).ln()/D(306)-D(119)*D(2).ln()/D(201960)-D(541)*D(P).ln()/D(201960)).exp()
 transport=(D(14000).ln()+D(12533)/D(459000)*(D(5)/D(4)).ln()).exp()
 admissible=((D(81127)*D(2).ln()+D(19068)*D(P).ln()-D(14616)*D(P+1).ln())/D(142783)).exp()/D(3)
 logratio38=D(201960)*(ordinary/D(8872)).ln()
 logratio40=D(459000)*(transport/D(14086)).ln()
 logratio41=D(142783)*(D(28885)/admissible).ln()
 gates=[{'equation':38,'scalar_value_display':str(ordinary),'scalar_cap':8872,'scalar_margin_display':str(D(8872)-ordinary),'normalized_ratio_display':str(ordinary/D(8872)),'raw_integer_LHS_over_RHS_display':str(logratio38.exp()),'exact_original_integer_gate_pass':21718**201960*2541935**660<8872**201960*2**119*P**541},
 {'equation':40,'threshold_display':str(transport),'chosen_threshold':14086,'threshold_margin_display':str(D(14086)-transport),'normalized_ratio_display':str(transport/D(14086)),'small_integer_side_over_large_display':str(logratio40.exp()),'exact_original_integer_gate_pass':14086**459000*4**12533>14000**459000*5**12533},
 {'equation':41,'admissible_C_display':str(admissible),'chosen_C':28885,'C_margin_display':str(admissible-D(28885)),'normalized_ratio_display':str(D(28885)/admissible),'raw_integer_LHS_over_RHS_display':str(logratio41.exp()),'exact_original_integer_gate_pass':(3*28885)**142783*(P+1)**14616<2**81127*P**19068}]
assert all(g['exact_original_integer_gate_pass'] for g in gates)
out={'status':'PASS_12_RATIONAL_MARGIN_BRACKETS_AND_3_EXACT_INTEGER_GATES','source_file':source.name,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'equations':gates,'scope':'The 12 endpoint margins are rigorous rational brackets of the saved upper enclosures, with fresh exact power comparisons. Equations 38/40/41 exact integer inequalities are freshly checked; their normalized margin decimals are display-only, not directed interval certificates. This is supplementary documentation and changes no paper file.'}
(O/'NUMERICAL_MARGINS.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
md=['# Numerical margins for the frozen paper','','The twelve rows below bound the expression in Lemma 9.1 using the stored rational upper bounds for U, W and the trace. Each interval was checked afresh by exact rational exponentiation; it is not obtained by rounding a displayed decimal.','','| Trace branch | Endpoint | Certified positive margin below 28885 |','|---|---|---|']
for r in rows:md.append(f"| {r['region']} | {r['point']} | [{r['margin_from_28885_enclosure'][0]}, {r['margin_from_28885_enclosure'][1]}] |")
md+=['','The tightest row is the count crossover x_m, with positive margin at least '+min(rows,key=lambda r:F(r['margin_from_28885_enclosure'][0]))['margin_lower_decimal']+'.','','## Scalar gates','','All three original integer inequalities were checked afresh with integer exponentiation. The following decimals are explanatory displays, not directed interval certificates.','','| Equation | Normalized interpretation | Displayed margin |','|---|---|---|',f"| (38) | 8872 minus the ordinary-cap scalar | {gates[0]['scalar_margin_display']} |",f"| (40) | 14086 minus 14000(5/4)^(12533/459000) | {gates[1]['threshold_margin_display']} |",f"| (41) | Largest admissible endpoint constant minus 28885 | {gates[2]['C_margin_display']} |",'','These margins quantify the current bound. They do not by themselves permit individual certified input constants to change: sensitivity requires rerunning the dependent comparisons.','']
(O/'NUMERICAL_MARGINS.md').write_text('\n'.join(md),encoding='utf-8')
print(json.dumps({'status':out['status'],'tightest_row_margin':min(rows,key=lambda r:F(r['margin_from_28885_enclosure'][0]))['margin_lower_decimal'],'scalar_gates':gates},indent=2),flush=True)
