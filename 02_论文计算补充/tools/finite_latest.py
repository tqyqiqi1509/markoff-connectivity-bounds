"""Numerical certificate checks for (10^82,10^84]; optional independent full CRT.

Uses the original arithmetic and whole-condition CRT functions. The current
mathematical rule index replaces historical model-admission and process logs.
"""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib,math,collections,itertools,bisect,argparse,time
HERE=Path(__file__).resolve().parent;PACKAGE=HERE.parent;BINDINGS={}
def remember(p):
 p=Path(p);h=hashlib.sha256(p.read_bytes()).hexdigest();BINDINGS[str(p.relative_to(ROOT))]=h;return h
def binding(p):return {'sha256':remember(p)}
def load(p):remember(p);return json.loads(Path(p).read_text())
def signed(n, d):
    return (n - 1) % (2 * d) == 0 or (n + 1) % (2 * d) == 0

def numbers(path):
    remember(path)
    result = [int(s) for s in Path(path).read_text().splitlines()]
    assert result == sorted(set(result)), str(path)
    return result

def seq_sha(seq):
    return hashlib.sha256(''.join(str(x) + '\n' for x in seq).encode()).hexdigest()

def boolean_check(plan):
    vertices = sorted(set(itertools.chain.from_iterable(e['d'] for e in plan['mixed_edges'])))
    assert vertices == plan['vertices'] and len(vertices) == 13
    bit = {d: 1 << i for i, d in enumerate(vertices)}
    clauses = [sum(bit[d] for d in e['d']) for e in plan['mixed_edges']]
    truth = [int(all(mask & clause for clause in clauses)) for mask in range(8192)]
    minimum = []
    for mask, yes in enumerate(truth):
        if yes and all(not truth[mask ^ (1 << j)] for j in range(13) if mask >> j & 1):
            minimum.append([d for j, d in enumerate(vertices) if mask >> j & 1])
    minimum.sort(key=lambda s: (len(s), s))
    assert minimum == plan['minimum'] and len(minimum) == 45
    minimal_masks = [sum(bit[d] for d in s) for s in minimum]
    coeff = []
    for mask in range(8192):
        assert bool(truth[mask]) == any(mask & m == m for m in minimal_masks)
        total, sub = 0, mask
        while True:
            total += truth[sub] * (-1 if (mask ^ sub).bit_count() % 2 else 1)
            if sub == 0:
                break
            sub = (sub - 1) & mask
        coeff.append(total)
    assert truth == plan['truth'] and coeff == plan['mobius_coefficients']
    nz = [(m, c) for m, c in enumerate(coeff) if c]
    assert len(nz) == 86
    assert [(b['mask'], b['weight']) for b in plan['branches']] == nz
    for index, b in enumerate(plan['branches']):
        assert b['index'] == index
        assert b['extra'] == [d for j, d in enumerate(vertices) if b['mask'] >> j & 1]
        assert b['minimal_index'] == (minimal_masks.index(b['mask']) if b['mask'] in minimal_masks else None)
    for mask in range(8192):
        assert sum(c for m, c in nz if mask & m == m) == truth[mask]
    return {'assignments': 8192, 'satisfying_assignments': sum(truth),
            'minimal_hitting_sets': 45, 'nonzero_Mobius_terms': 86,
            'truth_sha256': hashlib.sha256(bytes(truth)).hexdigest(),
            'coefficient_histogram': dict(sorted(collections.Counter(coeff).items()))}

def whole_blocks(ds):
    """Only small structural metadata; no reuse of either CRT enumerator."""
    groups = [{2 * d} for d in ds]
    changed = True
    while changed:
        changed = False
        for i in range(len(groups)):
            for j in range(i):
                if math.gcd(math.lcm(*groups[i]), math.lcm(*groups[j])) > 2:
                    groups[j].update(groups.pop(i))
                    changed = True
                    break
            if changed:
                break
    whole = [math.lcm(*s) for s in groups]
    assert sum(x % 4 == 0 for x in whole) == 1
    blocks = sorted(x if x % 4 == 0 else x // 2 for x in whole)
    assert all(math.gcd(a, b) == 1 for i, a in enumerate(blocks) for b in blocks[:i])
    assert all(b > 2 for b in blocks)
    return blocks

def sieve(limit):
    flags = bytearray([1]) * limit
    flags[0:2] = b'\0\0'
    for p in range(2, math.isqrt(limit - 1) + 1):
        if flags[p]:
            flags[p*p::p] = bytes(((limit - 1 - p*p) // p) + 1)
    return flags

def check_pocklington(db, flags):
    proved, count = set(), collections.Counter()
    for text in sorted(db, key=int):
        n = int(text)
        assert text == str(n) and n >= 2
        node = db[text]
        count[node['method']] += 1
        if node['method'] == 'trial_division':
            assert n < len(flags) and flags[n]
        else:
            assert node['method'] == 'pocklington_n_minus_1'
            F, seen = 1, set()
            for factor in node['factors']:
                q, e, a = int(factor['q']), factor['e'], factor['witness']
                assert q < n and q in proved and q not in seen
                assert type(e) is int and e >= 1 and type(a) is int
                seen.add(q)
                F *= q ** e
                assert (n - 1) % q**e == 0
                assert pow(a, n - 1, n) == 1
                assert math.gcd(n, pow(a, (n - 1) // q, n) - 1) == 1
            assert F == int(node['F']) and (n - 1) % F == 0 and F*F > n
        proved.add(n)
    return proved, dict(count)

def width(exponents):
    """Plain convolution, independent of the frozen engine's prefix method."""
    polynomial = [1]
    for e in exponents:
        product = [0] * (len(polynomial) + e)
        for degree, coefficient in enumerate(polynomial):
            for shift in range(e + 1):
                product[degree + shift] += coefficient
        polynomial = product
    assert polynomial == polynomial[::-1]
    assert max(polynomial) == polynomial[(len(polynomial) - 1) // 2]
    return max(polynomial)

def arithmetic_check(priority):
    source = ARITH / 'priority_scalar_certificates.json'
    rows = load(source)
    assert [int(r['p']) for r in rows] == priority
    window = load(ARITH / 'priority_window.json')
    assert window['input_sha256'] == binding(source)['sha256'] and window['L'] == 128
    assert window['records'] == [] and window['count'] == window['passed'] == 0 and window['remaining'] == []
    flags = sieve(1_000_001)
    db = load(ARITH / 'pocklington/certificate_dag.json')
    proved, dag_methods = check_pocklington(db, flags)
    small_primes = [q for q in range(2, 10000) if flags[q]]
    factors_used, report = set(), []

    def side_bounds(side, half):
        factors = [(int(q), e) for q, e in side['known_factors']]
        assert len(factors) == len(set(q for q, e in factors))
        R, B, K = int(side['remainder']), side['roughness_bound'], side['omega_tail_upper']
        assert B == 10000 and type(K) is int and K >= 0 and R >= 1
        product = 1
        for q, e in factors:
            assert type(e) is int and e > 0 and q >= 2
            assert (flags[q] if q < len(flags) else q in proved)
            factors_used.add(q)
            assert R % q
            product *= q ** e
        assert product * R == half and B ** K <= R < B ** (K + 1)
        assert all(R % q for q in small_primes)
        exponents = dict(factors)
        W = width(list(exponents.values()) + [1] * K)
        exponents[2] = exponents.get(2, 0) + 1
        Wo = width(list(exponents.values()) + [1] * K)
        S = Fraction(2 ** K)
        for q, e in factors:
            S *= 1 + e * Fraction(q - 1, q)
        A = Fraction(2 ** (len(exponents) + K)) * Fraction(B, B - 1) ** K
        for q in exponents:
            A *= Fraction(q, q - 1)
        return W, Wo, S, A, {'R': str(R), 'K': K, 'known_factor_count': len(factors)}

    for row in rows:
        p = int(row['p'])
        assert p > 10 ** 82 and p <= 10 ** 84 and p % 2
        minus, plus = side_bounds(row['minus'], (p - 1) // 2), side_bounds(row['plus'], (p + 1) // 2)
        M, Mo, S, A = minus[0] + plus[0], minus[1] + plus[1], minus[2] + plus[2], minus[3]
        safety = 9216 * A*A * Mo*Mo
        trace = S * M*M * (27*M*M + 6) * (27*M*M - 18*M - 3)
        capacity = M*M * (27*M**3 - 18*M*M - 3*M) ** 2
        assert row['m_upper'] == M and row['ordinary_m_upper'] == Mo
        for key, value in [('endgame_A_upper', A), ('trace_R_upper', S), ('safe4_rhs', safety), ('trace_rhs', trace)]:
            assert Fraction(*map(int, row[key])) == value
        assert int(row['capacity_rhs']) == capacity
        gates = {'safe4_pass': p > safety, 'trace_pass': 8*p > trace, 'capacity_pass': 8*p > capacity,
                 'conditional_prime_connected': p > safety and (8*p > trace or 8*p > capacity)}
        assert gates['conditional_prime_connected']
        assert all(type(row[k]) is bool and row[k] == v for k, v in gates.items())
        report.append({'p': str(p), 'gates': gates, 'm_upper': M, 'ordinary_m_upper': Mo,
                       'minus_partial': minus[4], 'plus_partial': plus[4],
                       'primality_claimed': False})
    large = sorted(q for q in factors_used if q >= len(flags))
    assert list(map(int, load(ARITH / 'ACCEPTED_LARGE_PRIMES.json'))) == large
    assert len(large) == 5 and len(proved) == 18104
    return {'rows':report,'fully_checked_dag_nodes':len(proved),'dag_methods':dag_methods}


def reconstruct(ds,lower,upper):
    moduli=[2*d for d in ds];groups=[]
    unseen=set(range(len(ds)))
    while unseen:
        component={unseen.pop()};front=list(component)
        while front:
            i=front.pop()
            neighbors={j for j in unseen if math.gcd(moduli[i],moduli[j])>2}
            unseen-=neighbors;component|=neighbors;front.extend(neighbors)
        groups.append(component)
    whole=[math.lcm(*(moduli[i] for i in c)) for c in groups]
    assert sum(m%4==0 for m in whole)==1
    blocks=sorted(m if m%4==0 else m//2 for m in whole)
    M=math.lcm(*moduli);assert math.prod(blocks)==M and 0<lower<upper<M
    assert all(math.gcd(a,b)==1 for i,a in enumerate(blocks) for b in blocks[:i])
    atoms=[M//b*pow(M//b,-1,b)%M for b in blocks]
    assert all(a%b==int(i==j) for i,a in enumerate(atoms) for j,b in enumerate(blocks))
    def half(selected):
        values=[0]
        for atom in selected:values=[(value+sign*atom)%M for value in values for sign in (-1,1)]
        assert len(values)==len(set(values));return sorted(values)
    L=half(atoms[::2]);R=half(atoms[1::2]);found=[]
    for a in L:
        for carry in (0,1):
            first=bisect.bisect_left(R,lower+1+carry*M-a)
            last=bisect.bisect_left(R,upper+1+carry*M-a)
            found.extend(a+b-carry*M for b in R[first:last])
    found.sort();assert found==sorted(set(found))
    assert all(lower<n<=upper and all(signed(n,d) for d in ds) for n in found)
    return found,{'modulus':str(M),'blocks':list(map(str,blocks)),'components':len(blocks),'left_size':len(L),'right_size':len(R)}
def main():
 global ROOT,DATA,ARITH
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reconstruct',action='store_true');args=ap.parse_args()
 assert not args.output.exists();start=time.monotonic();ROOT=args.project_root.resolve();DATA=ROOT/'computations/finite_extension_session/stage3_pool10_v2';ARITH=ROOT/'computations/finite_extension_session/stage3_pool10_arithmetic'
 plan=load(DATA/'PLAN.json');assert (int(plan['lower_exclusive']),int(plan['upper_inclusive']))==(10**82,10**84)
 rules=json.loads((PACKAGE/'indices/mixed_rules.json').read_text());bycase={x['case']:x for x in rules};generals=json.loads((PACKAGE/'indices/general_rank_59.json').read_text())
 assert sorted(plan['base'])==sorted(x['d'] for x in generals)
 for e in plan['mixed_edges']:
  k=bycase[e['case']];assert sorted(e['d'])==sorted(k['d']) and e['strict_p_gate']==k['gate']
 boolean=boolean_check(plan);values=numbers(DATA/'CANDIDATES.txt');assert len(values)==10369
 if args.reconstruct:
  union=set()
  for i,B in enumerate(plan['minimum']):
   v,meta=reconstruct(sorted(set(plan['base'])|set(B)),10**82,10**84);assert v==numbers(DATA/f'MINIMAL_{i:03d}_CANDIDATES.txt');union.update(v);print('CRT_BRANCH',i,len(v),flush=True)
  assert sorted(union)==values
 records=[json.loads(s) for s in (DATA/'CLASSIFICATION.jsonl').read_text().splitlines()];remember(DATA/'CLASSIFICATION.jsonl');assert [int(x['n']) for x in records]==values
 survivors=[];composite=0
 for n,row in zip(values,records):
  assert 10**82<n<=10**84 and n%2==1 and all(signed(n,d) for d in plan['base'])
  assert row['base_residues']==[n%(2*d) for d in plan['base']]
  masks=[sum(1<<j for j,d in enumerate(e['d']) if signed(n,d)) for e in plan['mixed_edges']];assert masks==row['mixed_clause_satisfied_masks'] and all(masks)
  power=((n-1)&-(n-1)).bit_length()-1;exponent=(n-1)>>power;chain=[pow(2,exponent<<j,n) for j in range(power)];passed=chain[0]==1 or n-1 in chain
  if passed:assert row['classification']=='MR_BASE2_SURVIVOR_NOT_PRIMALITY_CLAIM';survivors.append(n)
  else:
   assert row['classification']=='STRICT_COMPOSITE_BASE2' and row['witness']=={'base':2,'odd_exponent':str(exponent),'power_of_two':power,'residue_chain':list(map(str,chain))};composite+=1
 assert composite==9843 and survivors==numbers(DATA/'MR_SURVIVORS.txt') and len(survivors)==526
 mixed=[];direct=[]
 for n in survivors:
  hit=next((e for e in rules if n>e['gate'] and all(not signed(n,d) for d in e['d'])),None)
  (mixed if hit else direct).append(n)
 assert len(mixed)==520 and len(direct)==6 and direct==list(map(int,load(ARITH/'PRIORITY_CANDIDATES.json')))
 arithmetic=arithmetic_check(direct)
 out={'status':'PASS_LATEST_FINITE_NUMERICAL_CERTIFICATES','independent_full_CRT_reconstructed':args.reconstruct,'candidate_count':len(values),'composite_witnesses':composite,'mixed_applications':len(mixed),'direct_cases':len(direct),'boolean':boolean,'arithmetic':arithmetic,'candidate_sha256':seq_sha(values),'survivor_sha256':seq_sha(survivors),'input_sha256':BINDINGS,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Fresh integer, witness, rule-application and prime-DAG checks. Mixed and single-parameter rank determinants are separate exact certificate inputs.','seconds':time.monotonic()-start}
 args.output.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ['arithmetic','input_sha256','boolean']},indent=2))
if __name__=='__main__':main()
