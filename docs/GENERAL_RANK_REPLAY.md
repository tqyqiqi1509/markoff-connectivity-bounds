# General rank certificates: regeneration and verification

This guide supplements `indices/general_rank_inputs.json`. That index is the
authoritative list of all 57 parameters, the actual N, the strict characteristic
gate and the exact certificate path/hash. Do not substitute N=4d for a listed N.
All paths below are relative to a restored working root, not to `source/` in the
Git checkout. Use a separate writable working copy for expensive runs: the
original programs often write or resume fixed result names beside their inputs.
Python assertions must remain enabled. The original programs use Linux/WSL,
Python with python-flint 0.9.0, and for the C++ generators GMP and OpenMP.

## Selected maximal minors

For the original-v14 family the stored normalized matrix has complete rows and
row metadata recording all divided contents. The generic determinant checker is

```sh
python prior_restored/original_v14/markoff_project/agents/computational_mechanisms/verify_minor_certificates.py CERTIFICATE.json --flint
```

`CERTIFICATE.json` must be its actual indexed path. Without `--flint`, this entry
uses its separate pure-Python Bareiss determinant implementation. It recomputes
the selected determinants and checks content support and the rough gcd. It does
not regenerate the polynomial source. This command is not the d121/d125 pipeline.

The independent source generator for prime powers is `phi_stream.cpp` in that
same directory. Its actual interface and build command are:

```sh
g++ -O3 -std=c++17 -fopenmp phi_stream.cpp -lgmpxx -lgmp -o phi_stream
./phi_stream D new_raw.jsonl 2 N
python finalize_stream_matrix.py D new_raw.jsonl --prefix new --reference ORIGINAL_MATRIX.json --method minors
```

Here `D` and `N` are the indexed d,N. Run these commands in a disposable copy of
that source directory. The finalizer compares every regenerated normalized row
and content division with the original matrix when `--reference` is supplied.
The output raw file must be new for a regeneration: the C++ entry resumes rows
already present. For small cases `martin_kappa0.py D --nd N --prefix new` is the
direct polynomial implementation, also writing beside its own source file.
Composite d uses `phi_stream_composite.cpp`, `composite_parameters.py` and
`finalize_composite_matrix.py`; its filter is not the prime-power filter.

For d=22,34,38 use the local `verify_composite_certificate.py` in
`prior_restored/baragar_run/computations/modular_rank_campaign/even_composite_coupling/`;
for d=33,39 use its counterpart in `composite_coupling/`. The interface is
`python LOCAL_VERIFIER.py D`. The d46 and d58 source directories carry their
own composite certificates, filter audits and modular journals. For d58 the
specialized local command is:

```sh
python agents/attachment_audit/d58/verify_d58_certificate.py 58
```

These verifiers reconstruct normalization/filter data and validate determinant
residues under explicit Hadamard bounds. Saved journals are resumed; to request
new modular determinants, use a copied directory without the prior journal and
header. The snapshot in the repository/release must stay intact.

## Codimension-one block certificates

For the original-v14 prime block family, use the indexed wrapper with the generic
minor checker above; it dispatches to `verify_codimension_one_certificate.py`.
The latter can also be invoked directly on its nested certificate:

```sh
python prior_restored/original_v14/markoff_project/agents/computational_mechanisms/verify_codimension_one_certificate.py NESTED_CERTIFICATE.json --flint
```

This checks the primitive kernel, two determinants, all integer/rational linear
identities and residual numerator gcd. `codimension_one_certificate.py` is its
producer. A selected nonzero minor alone is not the whole certificate.

## d121 and d125: exact cofactor and first-minor proofs

These are separate formats, not ordinary prime codimension-one certificates.
Their indexed wrappers name the exact cofactor, frozen base-minor, solve and
independent-validation files. The d125 wrapper additionally names a classical
exact first-minor proof so the optional Arb determinant enclosure need not be a
mathematical premise. Both versions are retained with their original identities.

The implementation lives in the original-v14 project's
`agents/divisor_assisted_crt/`, `agents/resumable_exact_determinants/` and `code/`.
`produce_one_rhs.py` constructs and exactly checks the primitive Q,D equation;
`import_divisor_crt.py` reconstructs the base determinant from a certified divisor
and a sufficiently long modular stream. Its supported interface is:

```sh
python agents/divisor_assisted_crt/import_divisor_crt.py SOURCE_MATRIX.json --qd QD.json --crt-state CRT.json --output NEW_EXACT.json --recompute-residues --python-verify --export-base-minors NEW_BASE.json
```

Run that command from the restored original-v14 project root. It reports an
exact determinant only when its strict quotient reconstruction bound holds.
`partial_not_a_certificate` is not completion. The d125 classical input stream
is `agents/divisor_assisted_crt/d125_classical_exact.imported_crt.json`, and its
primitive equation input is `agents/divisor_assisted_crt/d125_QD.json`.

`code/one_cofactor_column_certificate.py` produces the replacement-minor gcd
certificate after the base determinant is available. The independent
`code/verify_one_cofactor_candidate.py` takes `--candidate`, `--source`, `--base`,
`--solve`, `--crt` and `--output`. For d121, every exact input is recorded under
those five keys in
`agents/resumable_exact_determinants/d121_cofactor_independent_validation.json`.
For d125, the final exact witness is
`results/d125_final_cofactor_integer_witness.json`, together with the classical
first-minor proof named by its indexed wrapper. No single universal command is
claimed to freshly regenerate both exceptional cases in this publication run.

## d131,137 and d139,149,151,157,163

Each local source directory contains `full_dD_matrix.json`, the universal block
certificate and a local `verify_universal_modular.py` accepting positional D.
For example:

```sh
python prior_restored/baragar_run/computations/modular_rank_campaign/verify_universal_modular.py 131
python new_uniform/d163/verify_universal_modular.py 163
```

The parallel wrapper exposed by the public tools is tested only for
139,149,151,157,163. Do not automatically apply that wrapper or its fixed source
hash to d131/d137 or the composite family. Its example for d163 is:

```sh
python agents/parallel_modular_verifier_round40/verify_universal_modular_parallel.py 163 --source-dir new_uniform/d163 --output-dir NEW_DIRECTORY --workers 4
```

It reconstructs each claimed determinant modulo proven auxiliary primes until
the modulus product exceeds twice the explicit Hadamard bound. Reading an old
completion record, checking file hashes, and redoing these determinants are
different checks.

## Verification scope

This packaging audit restored the 57 indexed roots and their direct scientific
path dependencies. It found and repaired missing original-v14 aliases and the
d121/d125 relative-path cofactor chains. It did not freshly regenerate all 57
source-polynomial families and all their determinants. The release's file hashes
establish exact artifact identity, not mathematical correctness by themselves.

## Parameter groups in this snapshot

- selected maximal minors: 4, 5, 7, 8, 9, 11, 13, 15, 16, 17, 19, 22, 23, 25, 27, 29, 31, 32, 33, 34, 35, 37, 38, 39, 41, 43, 46, 47, 49, 53, 58, 59, 61, 64, 81, 101.
- codimension one: 67, 71, 73, 79, 83, 89, 97, 103, 107, 109, 113, 127.
- special cofactor: 121, 125.
- bounded modular block: 131, 137, 139, 149, 151, 157, 163.
