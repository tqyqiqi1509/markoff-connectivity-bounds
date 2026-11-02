#!/usr/bin/env python3
"""Reconstruct det(A) from a primitive one-RHS divisor and old CRT residues.

No new modular determinant is generated except optional recomputation checks.
An unfinished quotient reconstruction is never an exact determinant certificate.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
from common import (atomic, crt, fmpz_mat, load_verified_qd, read_snapshot,
                    require, sha, stamp)


def reconstruct(source, qd_file, crt_file, output, recompute=False,
                pure_python=False, export_base=None):
    started = time.monotonic()
    protected = {Path(source).resolve(), Path(qd_file).resolve(), Path(crt_file).resolve()}
    require(Path(output).resolve() not in protected, "Output would overwrite an input")
    rows, binding, src, Q, D, qd, qd_sha, qd_verification = load_verified_qd(
        source, qd_file, pure_python)
    state, state_sha, state_payload = read_snapshot(crt_file)
    # Preserve the exact atomic snapshot read, even while the old producer advances.
    snapshot = Path(output).with_name(Path(output).stem + ".imported_crt.json")
    require(snapshot.resolve() not in protected, "Snapshot path would overwrite an input")
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    temporary = snapshot.with_name(snapshot.name + ".tmp")
    temporary.write_bytes(state_payload)
    temporary.replace(snapshot)
    reduced, c, reduction = crt.content_reduce(rows)
    legacy_binding = {k: binding[k] for k in ("source_file", "source_sha256", "shape", "row_indices")}
    require({k: v for k, v in state["binding"].items() if k != "source_file"}
            == {k: v for k, v in legacy_binding.items() if k != "source_file"},
            "Legacy CRT source/ordered-row binding mismatch")
    if c == 0:
        result = {"schema": "divisor_assisted_crt_v1", "status": "certified_exact_zero",
                  "determinant": "0", "binding": binding, "content_reduction": reduction,
                  "proof": "An exact zero row or column of the selected integer matrix",
                  "created_utc": stamp()}
        atomic(output, result)
        return result
    hadamard = crt.hadamard_data(reduced)
    require(int(state["content_factor"]) == c, "Legacy content factor mismatch")
    # Same validator as the original resumable engine: prime proofs, all records,
    # accumulated R/M, source identity, content factors and exact Hadamard data.
    old_r, old_m = crt.validate_state(
        state, legacy_binding, reduction, hadamard,
        recompute=recompute, A=fmpz_mat(reduced) if recompute else None)
    Dred = D // math.gcd(D, c)
    bound = (1 << hadamard["bound_bits"]) // Dred
    residue, modulus = 0, 1
    kept, skipped = [], []
    for record in state["records"]:
        p, dm = record["p"], record["det_mod_p"]
        if Dred % p == 0:
            require(dm == 0, "A skipped prime dividing Dred has an inconsistent nonzero residue")
            skipped.append({"p": p, "reason": "p divides Dred; division not invertible"})
            continue
        qm = dm * pow(Dred % p, -1, p) % p
        residue, modulus = crt.crt_step(residue, modulus, p, qm)
        kept.append({"p": p, "quotient_mod_p": qm})
    complete = modulus > 2 * bound
    result = {
        "schema": "divisor_assisted_crt_v1",
        "status": "certified_exact" if complete else "partial_not_a_certificate",
        "binding": binding, "QD_file": str(qd_file), "QD_sha256": qd_sha,
        "QD_verification": qd_verification,
        "old_crt_file": str(crt_file), "old_crt_snapshot_file": str(snapshot),
        "old_crt_snapshot_sha256": state_sha,
        "old_crt_status": state["status"],
        "old_crt_record_count": len(state["records"]),
        "modular_determinants_recomputed": recompute,
        "modular_provenance": "All retained modular determinants independently recomputed"
            if recompute else "Exact modular determinant records imported from the bound legacy run; not recomputed here",
        "content_factor": str(c), "content_reduction": reduction, "hadamard": hadamard,
        "primitive_denominator": str(D), "reduced_divisor": str(Dred),
        "quotient_absolute_bound": str(bound),
        "quotient_bound_bits": bound.bit_length(),
        "kept_records": kept, "skipped_records": skipped,
        "quotient_residue": str(residue), "quotient_modulus": str(modulus),
        "quotient_modulus_bits": modulus.bit_length(),
        "strict_termination_holds": complete,
        "termination": "product of kept distinct proved primes > 2 * floor(2^b / Dred)",
        "nonprobabilistic_reconstruction": True,
        "created_utc": stamp(), "elapsed_seconds": time.monotonic()-started,
    }
    if complete:
        quotient = residue if 2 * residue <= modulus else residue - modulus
        require(abs(quotient) <= bound, "Centered quotient exceeds its exact bound")
        reduced_det = Dred * quotient
        determinant = c * reduced_det
        require(reduced_det % old_m == old_r, "Recovered determinant contradicts imported CRT records")
        if state.get("status") == "certified_exact":
            require(int(state["determinant"]) == determinant, "Prior complete CRT disagrees")
        result.update(quotient=str(quotient), reduced_determinant=str(reduced_det),
                      determinant=str(determinant), determinant_nonzero=determinant != 0)
    require(sha(source) == binding["source_sha256"], "Source changed during import")
    require(sha(qd_file) == qd_sha, "QD checkpoint changed during import")
    require(sha(snapshot) == state_sha, "Saved imported CRT snapshot changed")
    atomic(output, result)
    if export_base is not None:
        require(complete, "Cannot export a base determinant from a partial reconstruction")
        require(Path(export_base).resolve() not in protected | {Path(output).resolve(), snapshot.resolve()},
                "Export would overwrite an input/proof file")
        require(not Path(export_base).exists(), "Use a new frozen base-certificate destination")
        atomic(export_base, {
            "status": "certified_first_minor_only",
            "matrix_sha256": binding["source_sha256"],
            "d": src.get("d") if isinstance(src, dict) else None,
            "minor_certificates": [{"row_indices": binding["row_indices"], "det": result["determinant"]}],
            "proof_file": str(output), "proof_sha256": sha(output),
            "provenance": "Deterministic divisor-assisted CRT; not a whole-rank or rough-gcd certificate",
        })
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("--qd", type=Path, required=True)
    ap.add_argument("--crt-state", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--recompute-residues", action="store_true")
    ap.add_argument("--python-verify", action="store_true")
    ap.add_argument("--export-base-minors", type=Path)
    a = ap.parse_args()
    r = reconstruct(a.source, a.qd, a.crt_state, a.output, a.recompute_residues,
                    a.python_verify, a.export_base_minors)
    print(json.dumps({k: r.get(k) for k in (
        "status", "old_crt_record_count", "quotient_bound_bits",
        "quotient_modulus_bits", "strict_termination_holds",
        "modular_determinants_recomputed", "elapsed_seconds")}), flush=True)


if __name__ == "__main__":
    main()
