#!/usr/bin/env python3
"""New exact scalar certificate for the sign-quotient / half-modulus interface.

All PASS decisions use integers and Fraction. Decimal only locates candidate
root brackets, each subsequently checked by exact integer powers. The large
positive and ordinary spectral integrals remain explicitly inherited inputs.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys
import time

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
SESSION = ROOT.parents[1]
BASE = ROOT.parent / "attachment_audit" / "extracted"
RESTORED = SESSION / "prior_restored"
BITS = 80
S = 1 << BITS


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upmul(a, b):
    return -(-(a * b) // S)


def frac(obj):
    return F(int(obj["numerator"]), int(obj["denominator"]))


def pack(a):
    return {"numerator": str(a.numerator), "denominator": str(a.denominator)}


def primes_before(n):
    return [q for q in range(2, n) if all(q % j for j in range(2, math.isqrt(q) + 1))]


def integer_root(n, k):
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            break
        x = y
    while (x + 1) ** k <= n:
        x += 1
    while x ** k > n:
        x -= 1
    assert x ** k <= n < (x + 1) ** k
    return x


def power_box(x, a, b):
    if a < 0:
        x, a = 1 / x, -a
    num = x.numerator ** a * S ** b
    den = x.denominator ** a
    lo = integer_root(num // den, b)
    assert lo ** b * den <= num < (lo + 1) ** b * den
    hi = lo + int(lo ** b * den != num)
    return F(lo, S), F(hi, S)


def prime_power_box(q, a, b):
    """q^(-a/b), exact dyadic bracket, using Decimal only as a locator."""
    with localcontext() as ctx:
        ctx.prec = 100
        lo = int(Decimal(S) * (Decimal(q) ** (-Decimal(a) / Decimal(b))))
    target = S ** b
    qa = q ** a
    assert qa * lo ** b <= target < qa * (lo + 1) ** b
    return lo, lo + 1


def trace_certificate():
    rows = []
    total = low = S
    ps = primes_before(401)
    for q in ps:
        wl, wu = prime_power_box(q, 829, 9180)
        sl, su = prime_power_box(q, 1, 18)
        # B(infinity) <= 1+(1-1/q)*s_upper/(1-s_upper).
        binf = S + -(-((q - 1) * su * S) // (q * (S - su)))
        e, wp, sp, sb, local = 0, S, S, 0, S
        values = [{"e": 0, "upper_integer": S}]
        while True:
            e += 1
            wp, sp = upmul(wp, wu), upmul(sp, su)
            sb += sp
            if wp * binf <= S * S:
                break
            be = S + -(-((q - 1) * sb) // q)
            value = upmul(wp, be)
            local = max(local, value)
            values.append({"e": e, "upper_integer": value})
            assert e < 1000
        total = upmul(total, local)
        if q >= 19:
            low = upmul(low, local)
        rows.append({"q": q, "w_lower_integer": wl, "w_upper_integer": wu,
                     "s_lower_integer": sl, "s_upper_integer": su,
                     "B_infinity_upper_integer": binf, "first_omitted_e": e,
                     "tail_w_power_upper_integer": wp, "local_upper_integer": local,
                     "finite_values": values})
    first = 401
    wl, wu = prime_power_box(first, 829, 9180)
    sl, su = prime_power_box(first, 1, 18)
    assert wu * (S + su) < S * S
    assert total < 159842 * S and low < 1826 * S
    return {"status": "PASS_FRESH_INTEGER_ROOTS_AND_ALL_EXPONENT_TAILS",
            "K_high_upper": pack(F(total, S)), "K_low_odd_upper": pack(F(low, S)),
            "caps": {"high": 159842, "low_odd": 1826, "joint_product": 159842},
            "rows": rows, "first_inactive_prime": first,
            "inactive_root_integers": {"w_lower": wl, "w_upper": wu, "s_lower": sl, "s_upper": su},
            "support_rule": "gcd((p-1)/2,(p+1)/2)=1; no shared factor 2"}


def negative_certificate():
    ps = [q for q in primes_before(73) if q >= 19]
    rows, total = [], F(1)
    for q in ps:
        rl, ru = prime_power_box(q, 1, 17)
        sl, su = prime_power_box(q, 1, 18)
        e = 1
        while q ** e * (q - 1) ** 306 < q ** 306:
            e += 1
        r, s = F(ru, S), F(sl, S)
        values = []
        for j in range(e):
            num = sum(r ** v for v in range(j + 1))
            den = 1 + F(q - 1, q) * sum(s ** v for v in range(1, j + 1))
            values.append(num / den)
        local = max(values)
        total *= local
        rows.append({"q": q, "first_omitted_e": e,
                     "r_lower_integer": rl, "r_upper_integer": ru,
                     "s_lower_integer": sl, "s_upper_integer": su,
                     "ratios_upper": [pack(x) for x in values], "local_upper": pack(local)})
    assert 73 ** 305 < 72 ** 306
    assert total < F(11, 10)
    return {"status": "PASS_FRESH_ODD_NEGATIVE_MOMENT", "constant": "11/10",
            "exact_product_upper": pack(total), "local_rows": rows,
            "prime_tail_from": 73, "scope": "all n with every prime divisor >=19; n=1 included",
            "no_division_of_upper_bounds": True}


def tau_certificate():
    rows, total = [], S
    for q in primes_before(449):
        # The adjacent ratio is decreasing in the exponent, so this locates
        # the exact maximum; no sampled or unproved exponent truncation.
        e = 0
        while (e + 2) ** 44 > (e + 1) ** 44 * q ** 5:
            e += 1
        wl, wu = prime_power_box(q, 5, 44)
        wp = S
        for _ in range(e):
            wp = upmul(wp, wu)
        local = (e + 1) * wp
        total = upmul(total, local)
        rows.append({"q": q, "maximizing_e": e, "w_lower_integer": wl,
                     "w_upper_integer": wu, "local_upper_integer": local})
    assert 449 ** 5 > 2 ** 44
    assert total < 2541933 * S
    return {"status": "PASS_FRESH_TAU_EULER_PRODUCT", "epsilon": "5/44",
            "single_integer_bound": 2541933, "coprime_pair_integer_bound": 2541934,
            "working_pair_bound": 2541935, "product_upper": pack(F(total, S)),
            "rows": rows, "first_inactive_prime": 449}


def main():
    started = time.monotonic()
    mapping = {x["original_path"]: x for x in json.loads((BASE / "PATH_MAP.json").read_text())}
    deps = {}

    def read_old(logical):
        r = mapping[logical]
        p = BASE / r["stored_path"]
        assert sha(p) == r["stored_sha256"] and p.stat().st_size == r["stored_bytes"]
        deps[logical] = {"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size}
        return json.loads(p.read_text())

    positive = read_old("latest/positive157/ACCEPTED_HPLUS.json")
    old63 = read_old("latest/root/bad_branch_transfer/old63branches_and_prerequisites.json")
    assert positive["integer_upper"] == 34195
    assert frac(positive["exact_producer_upper"]) < 34195
    assert positive["all_recomputed_tail_values"] == 39424980
    op = RESTORED / "baragar_current_run/weighted_window_certificate/weighted_159.json"
    oi = op.with_name("replay_weighted_159.json")
    hashes = {Path(x["path"]).name: x["sha256"] for x in old63["sources"]}
    for p in [op, oi]:
        assert sha(p) == hashes[p.name]
        deps[p.name] = {"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size}
    ordinary, replay = json.loads(op.read_text()), json.loads(oi.read_text())
    assert ordinary["status"] == "ARITHMETIC PASSED" and replay["status"] == "PASSED"
    rs = 1 << replay["fixed_bits"]
    prod_high = frac(ordinary["recursion"]["support_constants"][63]["C_all_two_upper"])
    prod_odd = frac(ordinary["recursion"]["support_constants"][0]["C_odd_upper"])
    replay_high = F(int(replay["recursion"]["support_records"][63]["C_all_scaled"]), rs)
    replay_odd_row = replay["recursion"]["support_records"][0]["rows"][-1]
    assert replay_odd_row["q"] == 2
    replay_odd = F(int(replay_odd_row["C_next_scaled"]), rs)
    assert prod_high < 21718 and replay_high < 21718
    assert prod_odd < 301 and replay_odd < 301
    for i in range(63):
        assert frac(ordinary["recursion"]["assignments"][i]["joint_upper"]) < 14556
        assert F(int(replay["recursion"]["assignment_upper_scaled"][i]), rs) < 14556

    trace = trace_certificate()
    negative = negative_certificate()
    tau = tau_certificate()
    P = 2 * 10 ** 155
    H, L, HP, LN = F(8750), F(301), F(34195), F(11, 10)
    KH, KL, PC = F(159842), F(1826), F(159842)
    C = 34300
    eps = delta = F(829, 9180)
    sigma, eta, bp = F(1, 17), F(1, 18), F(2, 15)
    kappa = (6 + 3 * eta) / (1 - 3 * sigma)
    alpha, a = eps * kappa + delta, 1 / kappa
    extra = eta + sigma * kappa
    alpha_up, extra_low = F(174, 227), F(41, 84)
    assert eps + bp * delta == F(829, 8100) and eps == delta
    assert kappa == F(629, 84) and alpha == F(591077, 771120)
    assert extra == F(125, 252) and extra > extra_low and alpha < alpha_up < 1
    assert 6 + 3 * extra == kappa
    # The original safe8 and cage gates remain conservative for half orders.
    rt = integer_root((P + 1) ** 5, 44)
    T = rt + int(rt ** 44 != (P + 1) ** 5)
    assert (T - 1) ** 44 < (P + 1) ** 5 <= T ** 44
    safe = 36864 * 2541933 ** 4 * 2541935 ** 2 * T ** 6
    cage = 256 * 2541933 ** 4 * T ** 4
    assert safe < P and cage < P
    # New transfer: d >=sqrt(2p)/(KD X^(5/44)), X=(p+1)/2.
    ds, de = sigma - F(1, 18), eps - F(4, 45)
    beta = F(5, 44) * ds - de
    assert ds == F(1, 306) and de == F(13, 9180)
    assert beta == F(-211, 201960)
    assert ds / 2 - beta == F(541, 201960)
    assert ds / 2 + beta == F(119, 201960)
    assert 21718 ** 201960 * 2541935 ** 660 < int(H) ** 201960 * 2 ** 119 * P ** 541
    # Complete compact-domain envelope, using outward algebraic root boxes.
    xl, xu = power_box(HP / H, 15, 2)
    y0, x1, x2 = L / LN, PC / KL, PC / (L / LN)
    assert 0 < x1 < x2 < xl <= xu < KH and y0 < KL
    _, hp_right = power_box(KH, -2, 15)
    pts = [("top_left", F(0), KL, H, L),
           ("top_hyperbola", x1, KL, H, L),
           ("low_crossover", x2, y0, H, L),
           ("high_crossover", xu, PC / xl, H, LN * PC / xl),
           ("right_edge", KH, PC / KH, HP * hp_right, LN * PC / KH)]
    vals = []
    for name, x, y, h, l in pts:
        _, z = power_box(x + y, 84, 629)
        v = (h + l) * z
        assert v < C
        vals.append({"point": name, "upper": pack(v), "display_only": float(v)})
    # All non63 support branches retain their old independent proof.
    p63 = 10 ** 154
    old_left = (729 * (2 * 181265 + 3)) ** 900 * 27 ** 420 * 14556 ** 6660 * (p63 + 1) ** 673
    assert old_left < (16 * p63) ** 900 and P >= p63
    # p < (3C)^kappa X^alpha /(8*2^extra), X=(p+1)/2.
    # Relax alpha upward and extra downward; then raise to19068.
    assert kappa * 19068 == 142783
    assert alpha_up * 19068 == 14616
    assert (3 + alpha_up + extra_low) * 19068 == 81127
    left = (3 * C) ** 142783 * (P + 1) ** 14616
    right = (1 << 81127) * P ** 19068
    assert left < right
    main_geometry = ROOT.parent / "literature_geometry_round2/EVEN_SPECTRAL_COVER_SIMPLIFICATION.md"
    assert main_geometry.is_file()
    geometry = []
    for p in [main_geometry, SESSION / "root_logs/SIGN_QUOTIENT_WRONSKIAN_ROOT.md",
              ROOT.parent / "literature_rank/SIGN_QUOTIENT_COUNTING_AUDIT.md"]:
        if p.exists():
            geometry.append({"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size})
    result = {"status": "PASS_NEW_HALF_MODULUS_SCALARS_WITH_EXPLICIT_GEOMETRIC_AND_SPECTRAL_INPUTS",
              "conditional_claim": "Under the accepted old Wronskian lemma and the even-spectral-coset assembly, every prime p>=2*10^155 has connected pure Vieta graph.",
              "P": str(P), "C": C, "H0": int(H), "L0_odd": int(L),
              "parameters": {k: str(v) for k, v in [("eps", eps), ("delta", delta), ("sigma", sigma), ("eta", eta),
                              ("positive_power", bp), ("kappa", kappa), ("alpha", alpha), ("extra_half_order", extra),
                              ("alpha_upper", alpha_up), ("extra_lower", extra_low)]},
              "trace": trace, "negative_odd": negative, "tau": tau,
              "inherited_ordinary_extractions": {"producer_high": pack(prod_high), "replay_high": pack(replay_high),
                                                 "producer_odd": pack(prod_odd), "replay_odd": pack(replay_odd)},
              "five_point_envelope": vals,
              "gates": {"original_safe8": True, "original_cage": True, "new_half_order_transfer": True,
                        "all_other63_old_supports_at_1e154": True, "endpoint": True},
              "prerequisites": {"T": str(T), "safe8_left": str(safe), "cage_left": str(cage)},
              "endpoint": {"left_formula": "(3*34300)^142783*(2*10^155+1)^14616",
                           "right_formula": "2^81127*(2*10^155)^19068",
                           "left_bits": left.bit_length(), "right_bits": right.bit_length()},
              "dependencies": deps, "geometry_inputs": geometry,
              "main_geometric_entry": str(main_geometry),
              "quotient_geometry_role": "Independent cross-check only; the direct even-spectral lift uses the already audited old Wronskian lemma.",
              "fresh_vs_inherited": "Trace Euler factors, odd negative moment, ordinary tau Euler factors, every scalar gate and continuous-domain envelope freshly recomputed; positive and ordinary spectral proofs/nodes inherited here, separate full pointwise rerun owned by root.",
              "source_sha256": sha(Path(__file__)), "seconds": time.monotonic() - started}
    (ROOT / "HALF_MODULUS_ENDPOINT_CERTIFICATE.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ["status", "conditional_claim", "P", "C", "H0", "gates", "endpoint", "seconds"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
