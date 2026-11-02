#!/usr/bin/env python3
"""Independent all-tail Arb/acb certificate for maximal truncated divisors.

No numerical integration or recurrence implementation from the root project is
imported.  Complex polynomials are summed term by term.  Every node of every
tail is serialized as an outward-rounded fixed-point integer.  Summation is
then exact integer arithmetic; only the certified Lipschitz correction uses
Arb.  Float values are restricted to elapsed-time diagnostics.
"""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from math import isqrt, lcm
from pathlib import Path
import struct
import sys
import time

import flint
from flint import arb, acb, ctx

sys.set_int_max_str_digits(0)


def primes_upto(n):
    flags = bytearray([1]) * (n + 1)
    flags[0:2] = b"\0\0"
    for q in range(2, isqrt(n) + 1):
        if flags[q]:
            flags[q*q:n+1:q] = b"\0" * ((n - q*q)//q + 1)
    return [i for i in range(2, n+1) if flags[i]]


def dyad(x):
    assert x.is_exact() and x.is_finite()
    m, e = x.man_exp()
    return {"mantissa": str(m), "exponent": int(e)}


def upper(x):
    return x.upper()


def ceil_upper(x):
    return int(x.upper().ceil().fmpz())


def invpow(q, x):
    return arb(q).root(x.denominator) ** (-x.numerator)


def finite_prime_limit(eps):
    k = (eps.denominator + eps.numerator - 1)//eps.numerator
    limit = 1 << k
    assert limit**eps.numerator >= 2**eps.denominator
    return limit


def initialize(eps, sig, panels):
    limit = finite_prime_limit(eps)
    active, omitted, metadata = [], [], []
    found_inactive = False
    for q in primes_upto(limit):
        w, r = invpow(q, eps), invpow(q, sig)
        assert 0 < w and w < 1 and 0 < r and r < 1
        if upper(w*(1+r)) <= 1:
            found_inactive = True
            omitted.append(q)
            continue
        assert not found_inactive, "active primes must be an initial segment"
        E, we = 1, w
        while not upper(we/(1-r)) < 1:
            E += 1
            we *= w
            assert E < 10000
        powers, rpower, sum0, sum1 = [arb(1)], arb(1), arb(1), arb(0)
        K, L = arb(1), arb(0)
        for e in range(1, E):
            powers.append(powers[-1]*w)
            rpower *= r
            sum0 += rpower
            sum1 += e*rpower
            K = max(K, upper(powers[-1]*sum0))
            L = max(L, upper(powers[-1]*sum1))
        active.append({"q": q, "w": w, "r": r, "wp": powers,
                       "K": K, "L": L, "first_discarded": E})
        metadata.append({"q": q, "first_discarded_exponent": E,
            "w_lower": dyad(w.lower()), "w_upper": dyad(w.upper()),
            "r_lower": dyad(r.lower()), "r_upper": dyad(r.upper()),
            "discarded_tail_upper": dyad(upper(we/(1-r))),
            "K_upper": dyad(K), "L_upper": dyad(L)})
    assert omitted
    # Monotonicity of x^-eps+x^-(eps+sig) handles all omitted real x,
    # independently of the broader integer cutoff above.
    q0 = omitted[0]
    inactive_test = upper(invpow(q0, eps)*(1+invpow(q0, sig)))
    assert inactive_test <= 1
    n = len(active)
    Ktails, corrections = [arb(1) for _ in range(n+1)], [arb(1) for _ in range(n+1)]
    Delta = arb.pi()/(2*panels)
    for j in range(n-1, -1, -1):
        Ktails[j] = upper(active[j]["K"]*Ktails[j+1])
        corrections[j] = upper((1+active[j]["L"]*Delta)*corrections[j+1])
    return active, Ktails, corrections, {
        "prime_limit": limit, "active_primes": metadata,
        "first_omitted_prime": q0, "first_omitted_test_upper": dyad(inactive_test),
        "omitted_primes_up_to_limit": omitted,
        "pi_lower": dyad(arb.pi().lower()), "pi_upper": dyad(arb.pi().upper())}


def integrate_all_tails(active, Ktails, corrections, panels, scale_bits, out, block_size):
    n, scale = len(active), 1 << scale_bits
    codec = struct.Struct("<" + "Q"*(n+1))
    assert Ktails[0]*scale < 2**64
    totals, blocksums, blocks = [0]*(n+1), [0]*(n+1), []
    raw_path = out.with_suffix(".nodes.u64.gz")
    digest = hashlib.sha256()
    pi, started = arb.pi(), time.monotonic()
    with gzip.open(raw_path, "wb", compresslevel=6) as stream:
        for node in range(panels):
            z = acb(0, pi*(2*node+1)/(2*panels)).exp()
            gs = []
            for row in active:
                rz = row["r"]*z
                term, poly, g = acb(1), acb(1), arb(1)
                for e in range(1, len(row["wp"])):
                    term *= rz
                    poly += term
                    g = max(g, upper(row["wp"][e]*poly.abs_upper()))
                gs.append(g)
            values, product = [scale]*(n+1), arb(1)
            for j in range(n-1, -1, -1):
                product *= gs[j]
                values[j] = ceil_upper(product*scale)
                assert scale <= values[j] < 2**64
            encoded = codec.pack(*values)
            stream.write(encoded)
            digest.update(encoded)
            for j, value in enumerate(values):
                totals[j] += value
                blocksums[j] += value
            if (node+1) % block_size == 0 or node+1 == panels:
                blocks.append({"first_node": node+1-min(block_size, (node%block_size)+1),
                               "last_node": node, "tail_sums": list(blocksums)})
                blocksums = [0]*(n+1)
                progress = {"stage": "integration", "completed_nodes": node+1,
                           "total_nodes": panels,
                           "seconds_display_only": round(time.monotonic()-started, 3)}
                print(json.dumps(progress), flush=True)
                out.with_suffix(".checkpoint.json").write_text(json.dumps(progress)+"\n")
    # Read back every saved node and check every integer total.
    readback = [0]*(n+1)
    read_digest = hashlib.sha256()
    with gzip.open(raw_path, "rb") as stream:
        for node in range(panels):
            data = stream.read(codec.size)
            assert len(data) == codec.size
            read_digest.update(data)
            values = codec.unpack(data)
            for j, v in enumerate(values):
                readback[j] += v
        assert stream.read(1) == b""
    assert readback == totals and read_digest.digest() == digest.digest()
    Wtails, tails = [], []
    for j in range(n+1):
        midpoint = arb(totals[j])/(panels*scale)
        bound = upper(midpoint*corrections[j])
        bound = min(bound, Ktails[j])
        Wtails.append(bound)
        tails.append({"q": active[j]["q"] if j < n else None,
                      "K_upper": dyad(Ktails[j]),
                      "correction_upper": dyad(corrections[j]),
                      "node_sum_integer": str(totals[j]),
                      "midpoint_display_only": midpoint.str(35),
                      "W_upper": dyad(bound), "W_display_only": bound.str(35)})
    return Wtails, {"panels": panels, "completed_nodes": panels,
        "node_scale_bits": scale_bits, "node_file": raw_path.name,
        "node_format": "gzip of little-endian uint64; rows are nodes, columns q-ascending tail indices including final empty tail",
        "row_bytes": codec.size, "uncompressed_bytes": panels*codec.size,
        "node_raw_sha256": digest.hexdigest(),
        "node_gzip_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
        "all_node_readback_verified": True, "tails": tails, "blocks": blocks,
        "seconds_display_only": time.monotonic()-started}


def recurse(active, Ktails, Wtails):
    n, Cnext, records = len(active), arb(1), []
    exact_two = None
    for j in range(n-1, -1, -1):
        row = active[j]
        q, w, r = row["q"], row["w"], row["r"]
        # This uses the very same certified Cnext in both sides, so it proves
        # a finite supremum for the monotone upper-bound recurrence.
        A, wa = 1, w
        while not upper(wa*(Ktails[j+1]+Cnext)) <= Cnext:
            A += 1
            wa *= w
            assert A < 10000
        omitted_upper = upper(wa*(Ktails[j+1]+Cnext))
        assert omitted_upper <= Cnext
        best, best_a, candidates = Cnext, 0, []
        wa, ra = arb(1), arb(1)
        qpow = 1
        for a in range(1, A):
            wa *= w
            ra *= r
            qpow *= q
            J, chosen, layers, E = Ktails[j+1], None, None, arb(1)
            for k in range(j+1, n):
                Q = active[k]["q"]
                h, Qpower = 1, Q
                while Qpower < qpow:
                    h += 1
                    Qpower *= Q
                assert Q**h >= qpow and (h == 1 or Q**(h-1) < qpow)
                trial = upper(E*h*Wtails[k])
                if trial < J:
                    J, chosen, layers = trial, Q, h
                # E for the next cut; never divide two outward upper bounds.
                E = upper(E*active[k]["K"])
            value = upper(wa*J + wa*ra*Cnext)
            candidates.append({"a": a, "cut_prime": chosen, "layers": layers,
                               "J_upper": dyad(J), "value_upper": dyad(value),
                               "value_display_only": value.str(35)})
            if q == 2 and a == 1:
                exact_two = value
            if value > best:
                best, best_a = value, a
        records.append({"q": q, "Cnext_upper": dyad(Cnext),
            "first_omitted_a": A, "omitted_upper": dyad(omitted_upper),
            "best_a": best_a, "C_upper": dyad(best),
            "C_display_only": best.str(35), "candidates": candidates})
        Cnext = best
    assert exact_two is not None
    joint = upper(Cnext+exact_two)
    return joint, {"C_all_upper": dyad(Cnext), "C_all_display_only": Cnext.str(40),
        "C_exact_v2_one_upper": dyad(exact_two),
        "C_exact_v2_one_display_only": exact_two.str(40),
        "C_joint_upper": dyad(joint), "C_joint_display_only": joint.str(40),
        "C_joint_integer_upper": ceil_upper(joint), "rows": records}


def euler_certificate(exponent, weighted):
    limit, value, local = finite_prime_limit(exponent), arb(1), []
    x, y = exponent.numerator, exponent.denominator
    for q in primes_upto(limit):
        d, h = (q, q-1) if weighted else (1, 1)
        e = 0
        while (d+(e+1)*h)**y > (d+e*h)**y*q**x:
            e += 1
        assert (d+(e+1)*h)**y <= (d+e*h)**y*q**x
        if e:
            assert (d+e*h)**y > (d+(e-1)*h)**y*q**x
        value *= (arb(d+e*h)/d)*invpow(q, exponent)**e
        local.append([q, e])
    return ceil_upper(value), {"exponent": str(exponent), "weighted": weighted,
        "prime_limit": limit, "local_exponents": local,
        "value_lower": dyad(value.lower()), "value_upper": dyad(value.upper()),
        "value_display_only": value.str(40)}


def endpoint(Cjoint, eps, sig, delta, power):
    KA, Adata = euler_certificate(delta, True)
    Ktau, taudata = euler_certificate(Fraction(5,44), False)
    assert KA <= 200496 and Ktau <= 2541933
    P, C = 10**power, ceil_upper(Cjoint)
    beta = 1-3*sig
    gamma27, gammaC, alpha = 6*sig/beta, Fraction(6)/beta, delta+6*eps/beta
    N = lcm(gamma27.denominator, gammaC.denominator, alpha.denominator)
    factor = 216*(Fraction(KA)+Fraction(3,2))
    lhs = factor.numerator**N * 27**int(N*gamma27) * C**int(N*gammaC) * (P+1)**int(N*alpha)
    rhs = factor.denominator**N * P**N
    T = ceil_upper(arb(P+1).root(44)**5)
    assert (T-1)**44 < (P+1)**5 <= T**44
    safe = 36864*Ktau**4*(Ktau+2)**2*T**6
    cage = 256*Ktau**4*T**4
    ratio = (arb(factor.numerator)/factor.denominator * arb(27)**(arb(gamma27.numerator)/gamma27.denominator)
             * arb(C)**(arb(gammaC.numerator)/gammaC.denominator)
             * arb(P+1)**(arb(alpha.numerator)/alpha.denominator)/P)
    checks = {"alpha_lt_one": alpha < 1, "exact_main": lhs < rhs,
              "Arb_main": ratio < 1, "safe8": safe < P, "primitive_cage": cage < P}
    return {"power": power, "C_joint_integer_upper": C, "K_A": KA, "K_tau": Ktau,
        "beta": str(beta), "alpha": str(alpha), "gamma27": str(gamma27),
        "gammaC": str(gammaC), "common_power": N,
        "main_left_sha256": hashlib.sha256(str(lhs).encode()).hexdigest(),
        "main_right_sha256": hashlib.sha256(str(rhs).encode()).hexdigest(),
        "main_ratio_upper": dyad(upper(ratio)), "main_ratio_display_only": ratio.str(40),
        "T": str(T), "safe8_left": str(safe), "primitive_cage_left": str(cage),
        "checks": checks, "euler_A": Adata, "euler_tau": taudata}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epsilon", default="4/45")
    ap.add_argument("--sigma", default="1/18")
    ap.add_argument("--delta", default="3/26")
    ap.add_argument("--panels", type=int, default=16384)
    ap.add_argument("--bits", type=int, default=192)
    ap.add_argument("--node-scale-bits", type=int, default=24)
    ap.add_argument("--block-size", type=int, default=512)
    ap.add_argument("--power", type=int, default=178)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    ctx.prec = args.bits
    eps, sig, delta = map(Fraction, (args.epsilon, args.sigma, args.delta))
    assert eps > 0 and 0 < sig < Fraction(1,3) and delta > 0
    assert args.panels > 0 and args.bits >= 96 and 0 <= args.node_scale_bits <= 32
    out = Path(args.output)
    out.parent.mkdir(exist_ok=True, parents=True)
    active, Ktails, corrections, parameters = initialize(eps,sig,args.panels)
    print(json.dumps({"stage":"setup", "active_primes":len(active),
        "candidate_exponents":sum(len(r["wp"])-1 for r in active),
        "K_all_display_only":Ktails[0].str(30),
        "correction_all_display_only":corrections[0].str(30)}),flush=True)
    Wtails, integration = integrate_all_tails(active,Ktails,corrections,args.panels,
                                              args.node_scale_bits,out,args.block_size)
    joint, recursion = recurse(active,Ktails,Wtails)
    target = endpoint(joint,eps,sig,delta,args.power)
    lower = endpoint(joint,eps,sig,delta,args.power-1)
    result = {"status":"PASSED" if all(target["checks"].values()) else "ENDPOINT FAILED",
        "epsilon":str(eps), "sigma":str(sig), "delta":str(delta),
        "bits":args.bits, "flint_version":flint.__version__,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters":parameters, "integration":integration, "recursion":recursion,
        "endpoint":target, "lower_endpoint_diagnostic":lower}
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"C_joint_integer_upper":recursion["C_joint_integer_upper"],
        "C_all_display_only":recursion["C_all_display_only"],
        "C_exact_v2_one_display_only":recursion["C_exact_v2_one_display_only"],
        "C_joint_display_only":recursion["C_joint_display_only"],
        "checks":target["checks"],"main_ratio_display_only":target["main_ratio_display_only"],
        "output":str(out)},indent=2),flush=True)
    if not all(target["checks"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
