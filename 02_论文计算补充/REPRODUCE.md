# Reproduce the computational steps

Use a working copy of this directory. Expanding data and rerunning tools creates files; fresh outputs containing timings need not have the same bytes as the included receipts.

Python 3.11 or later suffices for integer/Fraction checks. The Arb checks were actually run with Linux Python 3.14.4 and python-flint 0.9.0. The frequency generator uses four Linux fork workers; legacy source tools also use Linux resource limits. C++ source generators require their declared GMP/FLINT headers and libraries. Native binaries and runtimes are not bundled.

1. Expand inputs and restore a readable work tree:
```
python tools/prepare_data.py
python tools/materialize.py /path/to/work/input_tree
```
`INPUTS.json` is authoritative for exact versions. The materializer reports any historical path collision instead of silently overwriting an earlier version.

2. Recompute analytic inputs (Arb tools require python-flint):
```
python tools/euler_constants.py
python tools/q2_windows.py --output /path/to/work/Q2.json
python tools/fractional_moment.py --q2 /path/to/work/Q2.json --output /path/to/work/MOMENT.json
python tools/resource/verify_resource_rational.py
MARKOFF_WORKDIR=/path/to/work/input_tree python tools/verify_selected_moment_independent.py
python tools/frequency_envelope.py
python tools/kernel_cells.py
python tools/certify_local_fourth.py
python tools/signed_trace.py
python tools/endpoint.py --moment /path/to/work/MOMENT.json --trace tools/EDGEWORTH_TRACE_CERTIFICATE.json --output /path/to/work/ENDPOINT.json
```
The q2/q3 checks cover all exponents and their infinite tails. The frequency generator recomputes exactly the required 15104 cells on [0,32]. Signed trace checks all finite variance boxes and the infinite variance tails. The final endpoint tool checks 41 gates and twelve endpoints, without relying on historical model-admission labels.

3. Finite-range checks:
```
python tools/check_finite.py --project-root /path/to/work/input_tree --output /path/to/work/BASE.json
python tools/finite_latest.py --project-root /path/to/work/input_tree --output /path/to/work/LATEST.json --reconstruct
python tools/case101.py
g++ -O2 -std=c++17 tools/small_graphs.cpp -o /path/to/work/small_graphs
/path/to/work/small_graphs
```
The first command reconstructs 27985 integers and 26210 composite witnesses, parses full source rows, and checks the 90 saved unit-minor batches for 1775 candidates. It does not recompute their large determinants. The second reconstructs all 45 CRT branches through 10^84, checks 9843 composite witnesses, 520 mixed-rule applications and the six direct arithmetic cases with 18104 prime-DAG nodes. Case101 freshly reconstructs its full source matrix and computes two 184-by-184 integer determinants. The small graph program checks all 211 primes from 5 through 1304.

For the preceding finite-extension intervals and the full general/mixed determinant replays, consult THEOREM_INPUTS.json and the exact original source/checker paths restored from INPUTS.json. Their integer sources, row normalization, selected minors, prime moduli, reconstruction bounds and journals are provided. Legacy command-line paths may need rebasing; retain the original bytes and record a path-only adapter. The sources are supplied to permit the mathematical checks, not to recreate discarded process-management state.
