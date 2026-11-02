#!/usr/bin/env python3
"""One fixed endpoint, independent Python integer/Fraction arithmetic.

No producer module, Arb, FLINT, floating-point log/root, LP cap, or q cap.
Every decision is exact. Decimal diagnostics are output only.
This is a conditional mathematical replay, not root acceptance of a new P.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback

sys.set_int_max_str_digits(0)

P = 8 * 10**152
C, H = 28885, Q(8872)
K, KL, LOW, CC = Q(159842), Q(662), Q(135), Q(26, 25)
CM, CE = Q(34195), Q(28924)
DELTA = Q(829, 9180) - Q(63, 1000)
LAMBDA, NU = Q(43, 1000), Q(61, 10000)
A, B = Q(9391, 610), Q(11200, 61)
ROOT_BITS, LOG_BITS, LOG_TERMS = 128, 192, 48

REQUIRED = {
    "control/RELAXED_BUDGET_TRACE_ROOT_ACCEPTANCE.json": (
        "87325c562c145910d00b6a490f7994c1c5dddcfd8839d7e151df52d32b6fc113",
        "ACCEPTED_GENERAL_RELAXED_BUDGET_TRACE_NO_NEW_P_ENDPOINT"),
    "control/RESOURCE_GAUSSIAN_SIMPLIFICATION_ROOT_ACCEPTANCE.json": (
        "f9f37783b2e0153266548bfb726510f1413fccc23efaae7694a1c37fc65141fe",
        "ACCEPTED_SIMPLIFIED_INTERNAL_PROOF_CHAIN_SAME_P_ENDPOINT"),
    "control/SUPPORT_SIMPLIFICATION_ROOT_ACCEPTANCE.json": (
        "bb177bf8e52ecfde06e3ab7a85f7dc77f7b01bda540669d26d6e78fde0ae473c",
        "ACCEPTED_SUPPORT_SIMPLIFICATION_NO_NEW_P_ENDPOINT"),
}

LOG_WITNESSES = []
ROOT_WITNESSES = []


def fraction(q):
    q = Q(q)
    return {"numerator": str(q.numerator), "denominator": str(q.denominator)}


def file_sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def integer_digest(n):
    assert isinstance(n, int) and n >= 0
    raw = n.to_bytes((n.bit_length() + 7) // 8 or 1, "big")
    return {"bits": n.bit_length(), "sha256_unsigned_big_endian": hashlib.sha256(raw).hexdigest()}


def grid_floor(q, bits):
    q = Q(q)
    scale = 1 << bits
    return Q((q.numerator * scale) // q.denominator, scale)


def grid_ceil(q, bits):
    return -grid_floor(-Q(q), bits)


def exact_floor_root(n, degree):
    """Integer Newton iteration followed by an explicit exact root bracket."""
    assert n >= 0 and degree > 0
    if n == 0:
        return 0
    if degree == 1:
        return n
    if degree == 2:
        root = math.isqrt(n)
    else:
        root = 1 << ((n.bit_length() + degree - 1) // degree)
        while True:
            next_root = ((degree - 1) * root + n // root**(degree - 1)) // degree
            if next_root >= root:
                break
            root = next_root
        while root**degree > n:
            root -= 1
        while (root + 1)**degree <= n:
            root += 1
    assert root**degree <= n < (root + 1)**degree
    return root


@lru_cache(maxsize=None)
def power_bounds(q, numerator, denominator):
    q = Q(q)
    assert q > 0 and denominator > 0
    target = q**numerator
    scale = 1 << ROOT_BITS
    radicand = (target.numerator * scale**denominator) // target.denominator
    integer_low = exact_floor_root(radicand, denominator)
    lower, upper = Q(integer_low, scale), Q(integer_low + 1, scale)
    assert lower**denominator <= target < upper**denominator
    ROOT_WITNESSES.append({
        "input": fraction(q), "power_numerator": numerator,
        "power_denominator": denominator, "lower": fraction(lower), "upper": fraction(upper),
        "bracket_checked_exactly": True,
    })
    return lower, upper


def atanh_log_bounds(reduced):
    """log(reduced), 1 <= reduced <= 2, positive 48-term sum and tail."""
    reduced = Q(reduced)
    assert 1 <= reduced <= 2
    z = (reduced - 1) / (reduced + 1)
    zz = z * z
    term = z
    partial = Q(0)
    for j in range(LOG_TERMS):
        partial += term / (2 * j + 1)
        term *= zz
    lower = 2 * partial
    tail = 2 * term / ((2 * LOG_TERMS + 1) * (1 - zz))
    assert tail >= 0
    return lower, lower + tail


LOG_TWO = atanh_log_bounds(Q(2))


@lru_cache(maxsize=None)
def log_bounds(q):
    q = Q(q)
    assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    reduced = q / (Q(2)**k)
    while reduced < 1:
        reduced *= 2
        k -= 1
    while reduced >= 2:
        reduced /= 2
        k += 1
    assert 1 <= reduced < 2 and q == Q(2)**k * reduced
    rlow, rhigh = atanh_log_bounds(reduced)
    if k >= 0:
        low, high = rlow + k * LOG_TWO[0], rhigh + k * LOG_TWO[1]
    else:
        low, high = rlow + k * LOG_TWO[1], rhigh + k * LOG_TWO[0]
    lower, upper = grid_floor(low, LOG_BITS), grid_ceil(high, LOG_BITS)
    assert lower <= low <= high <= upper
    LOG_WITNESSES.append({
        "input": fraction(q), "power_of_two": k, "reduced": fraction(reduced),
        "terms": LOG_TERMS, "rounding_bits": LOG_BITS,
        "lower": fraction(lower), "upper": fraction(upper),
        "tail_formula": "2*z^(2*N+1)/((2*N+1)*(1-z^2)), z=(r-1)/(r+1)",
    })
    return lower, upper


def decimal_display(q):
    q = Q(q)
    with localcontext() as ctx:
        ctx.prec = 32
        return str(Decimal(q.numerator) / Decimal(q.denominator))


def endpoint_display(u, w, trace):
    with localcontext() as ctx:
        ctx.prec = 32
        def dq(x):
            return Decimal(x.numerator) / Decimal(x.denominator)
        value = (Decimal(1510) * dq(u).ln() + Decimal(377) * dq(w).ln()
                 + Decimal(252) * dq(trace).ln()) / Decimal(1887)
        return str(value.exp())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--project-root", type=Path,
                        default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / "FRACTION_CONTINUOUS_BUDGET_REVIEW.json"
    if target.exists():
        raise SystemExit("Refusing to overwrite an existing review; choose a new --output-dir")
    start = time.monotonic()
    result = {
        "status": "RUNNING", "root_acceptance_claimed": False,
        "conjecture": "OPEN", "P": str(P), "C": C, "H": str(H),
        "input_parameters": {"K": str(K), "KL": str(KL), "L": str(LOW),
                             "c": str(CC), "Cm": str(CM), "CE": str(CE),
                             "delta_minus_epsilon": str(DELTA),
                             "cuts": [128, 10000, 14086]},
        "arithmetic": {"decision_backend": "Python integers and fractions.Fraction only",
                       "root_bits": ROOT_BITS, "log_bits": LOG_BITS, "log_terms": LOG_TERMS,
                       "decimal_fields": "display_only, never used for decisions",
                       "producer_imports": False, "round49_access": False},
        "dependencies": [], "gates": [], "endpoint_rows": [],
        "scope": "Fixed endpoint replay conditional on bound inherited graph/Euler/BE/mass/probability premises. Root alone may accept new p range.",
        "script_sha256": file_sha(__file__),
    }

    def persist():
        result["log_witnesses"] = LOG_WITNESSES
        result["root_witnesses"] = ROOT_WITNESSES
        result["elapsed_seconds_display_only"] = time.monotonic() - start
        target.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")

    def check(label, condition, **details):
        assert isinstance(condition, bool)
        result["gates"].append({"label": label, "pass": condition, **details})
        if not condition:
            result["status"] = "FAILED_FIXED_CANDIDATE_STOPPED"
            persist()
            raise AssertionError(label)

    try:
        for relative, (digest, status) in REQUIRED.items():
            path = args.project_root / relative
            actual = file_sha(path)
            check("required_source_hash:" + relative, actual == digest,
                  expected=digest, observed=actual)
            document = json.loads(path.read_text())
            check("required_source_status:" + relative, document["status"] == status)
            result["dependencies"].append({"path": relative, "sha256": actual,
                                           "status": document["status"]})
            for child, expected in sorted(document.get("bindings", {}).items()):
                observed = file_sha(args.project_root / child)
                check("inherited_binding:" + child, observed == expected,
                      expected=expected, observed=observed)

        check("exact_budget_delta", DELTA == Q(12533, 459000) and 0 < DELTA < Q(1, 4))
        check("exact_density_coefficient", A == Q(15, 2) + Q(28, 25) * LAMBDA / NU)
        check("exact_budget_coefficient", B == Q(28, 25) / NU)
        check("gaussian_density_rational_replacement", Q(18)**2 < Q(15, 2)**2 * 6)
        profile_constant = 4 * A**3 / (27 * B)
        check("whole_variance_profile_constant_lt_3", profile_constant < 3,
              value=fraction(profile_constant), display_only=decimal_display(profile_constant))
        check("middle_budget_transport", Q(4, 5) > Q(14, 15)**4)
        high_lhs = 14086**459000 * 4**12533
        high_rhs = 14000**459000 * 5**12533
        check("high_budget_transport", high_lhs > high_rhs,
              left=integer_digest(high_lhs), right=integer_digest(high_rhs))

        support = json.loads((args.project_root / "control/SUPPORT_SIMPLIFICATION_ROOT_ACCEPTANCE.json").read_text())
        check("support_scope_and_caps", P > 696
              and support["new_support_modulus"] == 12939386460
              and support["low_half_caps"] == {"B": 662, "count": 135, "negative_ratio": "26/25"})
        mods = [44, 66, 30, 70, 78, 68, 76, 92, 116]
        cumulative = mods[0]
        for modulus in mods[1:]:
            check("support_sign_coupling:" + str(modulus), math.gcd(cumulative, modulus) > 2)
            cumulative = math.lcm(cumulative, modulus)
        check("support_lcm", cumulative == 12939386460)

        ordinary_lhs = 21718**201960 * 2541935**660
        ordinary_rhs = int(H)**201960 * 2**119 * P**541
        check("ordinary_high_cap", ordinary_lhs < ordinary_rhs,
              left=integer_digest(ordinary_lhs), right=integer_digest(ordinary_rhs))
        radicand = (P + 1)**5
        floor_t = exact_floor_root(radicand, 44)
        ceiling_t = floor_t + (floor_t**44 < radicand)
        check("ceiling_bracket", (ceiling_t - 1)**44 < radicand <= ceiling_t**44,
              T=str(ceiling_t))
        safe_lhs = 36864 * 2541933**4 * 2541935**2 * (ceiling_t + 1)**6
        cage_lhs = 256 * 2541933**4 * (ceiling_t + 1)**4
        check("safe8_with_T_plus_1", safe_lhs < P, left=str(safe_lhs), right=str(P))
        check("cage_with_T_plus_1", cage_lhs < P, left=str(cage_lhs), right=str(P))
        check("ceil_uniform_safe_power", 6 * Q(5, 44) - 1 == -Q(7, 22))
        check("ceil_uniform_cage_power", 4 * Q(5, 44) - 1 == -Q(6, 11))
        check("final_ratio_uniform_derivative", 14616 < 19068)

        e_lower, e_upper = power_bounds(CE / H, 15, 2)
        m_lower, m_upper = power_bounds(CM / H, 15, 2)
        xh = K / KL
        xl = K * CC / LOW
        check("complete_crossover_order", 0 < 128 < xh < xl < e_lower <= e_upper
              < 10000 < 14086 < m_lower <= m_upper < K)
        result["crossovers"] = {"horizontal_to_hyperbola": fraction(xh),
                                "low_cap": fraction(xl),
                                "mass_lower": fraction(e_lower), "mass_upper": fraction(e_upper),
                                "count_lower": fraction(m_lower), "count_upper": fraction(m_upper)}
        logp_lower, logp_upper = log_bounds(Q(P, 2))
        logc_lower, logc_upper = log_bounds(Q(43540))
        d128_lower = log_bounds(Q(128))[0] + DELTA * logp_lower - logc_upper
        check("D128_strictly_gt_3", d128_lower > 3,
              D_lower=fraction(d128_lower), display_only=decimal_display(d128_lower))
        mix_middle = Q(2, 5) + Q(3, 5) * K / (10000**2 + K)
        mix_high = Q(177, 500) + (1 - Q(177, 500)) * K / (14086**2 + K)
        check("mixed_middle_total_cap", mix_middle < Q(41, 100), value=fraction(mix_middle))
        check("mixed_high_total_cap", mix_high < Q(355, 1000), value=fraction(mix_high))

        points = {
            "zero_limit": (Q(0), Q(0)), "cut_128": (Q(128), Q(128)),
            "horizontal_to_hyperbola": (xh, xh), "low_cap": (xl, xl),
            "mass_crossover": (e_lower, e_upper), "cut_10000": (Q(10000), Q(10000)),
            "cut_14086": (Q(14086), Q(14086)),
            "count_crossover": (m_lower, m_upper), "right_edge": (K, K),
        }
        regions = [
            ("trivial", ["zero_limit", "cut_128"], None),
            ("continuous_profile", ["cut_128", "horizontal_to_hyperbola", "low_cap", "mass_crossover", "cut_10000"], None),
            ("middle_41_100", ["cut_10000", "cut_14086"], Q(41, 100)),
            ("high_355_1000", ["cut_14086", "count_crossover", "right_edge"], Q(355, 1000)),
        ]
        result["continuous_domain_certificate"] = {
            "domain": "0<x<=K, 0<y<=KL, x*y<=K",
            "y_boundary": "min(KL,K/x)",
            "method": "monotonicity in y, then convexity of log F in log x on every cap branch",
            "profile_log_second_derivative": "1/(2*(log(x)+c0)^2)>0",
            "trace_jumps_both_sides": [128, 10000, 14086],
            "expected_endpoint_rows": 12,
            "zero_limit_included": True,
        }

        for region, labels, cap in regions:
            for label in labels:
                xlow, xupper = points[label]
                row = {"region": region, "point": label,
                       "x_lower": fraction(xlow), "x_upper": fraction(xupper)}
                if xlow == 0:
                    yupper = KL
                    uupper = wupper = H + LOW
                    traceupper = KL
                else:
                    yupper = min(KL, K / xlow)
                    power_upper = power_bounds(xlow, -2, 15)[1]
                    lowpart = min(LOW, CC * yupper)
                    uupper = min(H, CM * power_upper) + lowpart
                    wupper = min(H, CE * power_upper) + lowpart
                    if region == "continuous_profile":
                        dlow = log_bounds(xlow)[0] + DELTA * logp_lower - logc_upper
                        check("profile_D_gt_3:" + label, dlow > 3, D_lower=fraction(dlow))
                        gamma_upper = power_bounds(3 / dlow, 1, 2)[1]
                        check("profile_sqrt_strictly_lt_1:" + label, gamma_upper < 1)
                        traceupper = xupper * gamma_upper + yupper
                        row["D_lower"] = fraction(dlow)
                        row["sqrt_3_over_D_upper"] = fraction(gamma_upper)
                    elif cap is None:
                        traceupper = xupper + yupper
                    else:
                        traceupper = cap * (xupper + yupper)
                left = uupper.numerator**1510 * wupper.numerator**377 * traceupper.numerator**252
                right = C**1887 * uupper.denominator**1510 * wupper.denominator**377 * traceupper.denominator**252
                passed = left < right
                row.update({"y_upper": fraction(yupper), "U_upper": fraction(uupper),
                            "W_upper": fraction(wupper), "trace_upper": fraction(traceupper),
                            "total_constant_cap": None if cap is None else fraction(cap),
                            "integer_left": integer_digest(left), "integer_right": integer_digest(right),
                            "exact_gate_pass": passed,
                            "F_display_only": endpoint_display(uupper, wupper, traceupper)})
                result["endpoint_rows"].append(row)
                check("geometric_endpoint:" + region + ":" + label, passed)
        check("complete_endpoint_row_count", len(result["endpoint_rows"]) == 12)

        final_left = (3 * C)**142783 * (P + 1)**14616
        final_right = 2**81127 * P**19068
        check("final_endpoint_integer_gate", final_left < final_right,
              left=integer_digest(final_left), right=integer_digest(final_right))
        result["status"] = "PASS_FIXED_CONTINUOUS_BUDGET_ENDPOINT_PURE_INTEGER_FRACTION"
        result["new_probability_method_claimed"] = False
        result["removed_mathematical_dependencies"] = ["old_LP_effective_cap", "old_q_cap", "old_LP_filling"]
        result["retained_mathematical_dependencies"] = [
            "nine_rank_filters_and_support_bridge", "original_graph_and_Wronskian_chain",
            "actual_same_integer_Euler_budget_43540", "non_iid_BE_constant_14_25",
            "actual_mass_and_count_bounds", "accepted_general_2_5_and_177_500_trace_bounds"]
        persist()
        print(json.dumps({"status": result["status"], "P": str(P), "C": C,
                          "gates": len(result["gates"]), "endpoint_rows": len(result["endpoint_rows"]),
                          "output": str(target), "script_sha256": result["script_sha256"]}, indent=2), flush=True)
        for row in result["endpoint_rows"]:
            print(row["region"], row["point"], row["F_display_only"], row["exact_gate_pass"], flush=True)
        return 0
    except Exception as error:
        if result["status"] == "RUNNING":
            result["status"] = "ERROR_STOPPED_NO_ACCEPTANCE"
        result["exception"] = {"type": type(error).__name__, "message": str(error),
                               "traceback": traceback.format_exc()}
        persist()
        print(json.dumps({"status": result["status"], "error": str(error), "output": str(target)}, indent=2), flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
