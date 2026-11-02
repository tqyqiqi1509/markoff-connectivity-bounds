"""Exact additional gates and exponent identities used by the article."""
from fractions import Fraction as Q
from pathlib import Path
import json

P=8*10**152
T=237329177916957586
assert (T-1)**44 < (P+1)**5 <= T**44
assert 256*2541933**4*(T+1)**4 < P
assert 36864*2541933**4*(2*2541933)**2*(T+1)**6 < P
s,e,delta=Q(1,17),Q(1,18),Q(829,9180)
k=(9*s+1+2*e)/(1-3*s)
assert (8+2*k)/2==Q(1510,252)
assert (1+k)/2==Q(377,252)
assert (11+3*k)*delta/2==Q(591077,771120)
assert (9+3*k)/2==Q(629,84)
assert (5+k)/2==3+Q(125,252)
assert Q(5,44*306)+Q(4,45)-delta==-Q(211,201960)
assert Q(591077,771120)<Q(174,227)
assert Q(125,252)>Q(41,84)
result={'status':'PASS','P':str(P),'ceiling':T,'coarser_safe_count':'2 * 2541933',
        'scope':'Exact scalar gates and symbolic rational exponent identities only.'}
print(json.dumps(result,indent=2))
