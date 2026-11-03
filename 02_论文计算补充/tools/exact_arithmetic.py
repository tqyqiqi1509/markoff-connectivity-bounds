"""Exact interval arithmetic reused from the accepted round50 checker.
Original source SHA256: 87f39a975e547d03 (full SHA in INPUT_BINDINGS.json).
This module retains that implementation; it is not an independent algorithm.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib, math, sys
sys.set_int_max_str_digits(0)
ROOT_BITS, LOG_BITS, LOG_TERMS = 128, 192, 48
LOG_WITNESSES=[]
ROOT_WITNESSES=[]
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


