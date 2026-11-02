# Independent portability review

Review date: 2026-11-02. Scope: public-supplement selection and build scripts, original data dependencies, and the finite-check adapter. The manuscript was read only and was not changed. No repository was created and no files were uploaded by this reviewer.

## Findings repaired in the builder

1. **Wrong root for finite source inspection.** The former adapter changed only the global `ROOT`. The archived `inspect_source` function also reads global `OLD` and `NEW`, initialized from the reviewer's location. Moving that reviewer to `tools/` therefore made it look in the wrong directories. The adapter now updates `ROOT`, `OLD`, and `NEW` together before any functions are called.
2. **Two source-ready records require historical generator ELF identities.** Omitting only the `independent_modular_det` binary in each batch was insufficient. `d167_SOURCE_READY.json` and `d173_SOURCE_READY.json` also bind `new_pointwise/d*/phi_stream_partitioned`. These are compiler-produced executable bytes, excluded from the public source package. The adapter now recognizes only these two executable names, checks the corresponding original C++ source against its fixed digest, and writes each omitted executable identity to the result. It does not fabricate a digest, claim a matching binary rebuild, or claim a determinant was recomputed.
3. **The omission must be visible in the output.** The result now lists every omitted historical executable identity, its old digest, and the original source file/digest that was actually checked. All source-row, candidate, compositeness-witness, congruence, unit-determinant-receipt, and C++ source checks remain.

These changes are in `outputs/bargar_build_public.py`. The original archived scripts and certificates were not edited.

## Actual checks performed

The finite adapter was generated from the revised build-script fragment and executed against the original read-only input root, with a new output file outside that root.

- Complete CRT reconstruction: 27,985 candidates in `(1304,10^66]`.
- Compositeness witnesses recomputed: 26,210.
- Pointwise source-bound saved records checked: 1,775, in 90 batches.
- Historical executable identities explicitly omitted: 92 (90 batch determinant executables and 2 source generators).
- Result: PASS, about 20.9 seconds on this host.

The exact result is `PORTABLE_FINITE_FRESH_TEST.json`. Its scope is saved-certificate checking with source-row parsing and fresh CRT/compositeness arithmetic. It is **not** a fresh computation of the 1,775 determinants or the 57 general rank certificates.

Negative tests in `PORTABILITY_NEGATIVE_TESTS.json` confirm that changing a bound C++ source is rejected, and that a missing executable with an unrecognized name is not silently exempted.

## Selection issues reported to the coverage reviewer

The selected `certify_selected_moment.py` and `certify_size_trace.py` dynamically import `agents/analytic_window_round4/certify_trace_concentration.py`. That file was absent from the inspected selection. This is a genuine source-closure issue: either include it and its helper dependencies, or omit the old producers if the supported independent replay completely replaces their role. The coverage reviewer is handling this selection repair.

`prelaunch_history/finite_campaign.py`, `probe_gaussian_fine_variance.py`, and `probe_gaussian_tail_refinement.py` were also selected. Their names and interfaces indicate historical or exploratory work; they require an explicit proof-use justification or removal. A dependency map should distinguish a program actually required to regenerate a certificate from a historical path merely recorded inside that certificate.

The selector parses only JSON files below 3 MB and uses string/path heuristics. This is a useful discovery method but does not prove dependency closure. Large matrix/certificate JSON and dynamically constructed filenames must be checked separately. Missing references must be reported, not silently treated as irrelevant.

## Required release gates

Before calling the supplement portable or complete:

1. Finish all release ZIPs, produce `PATH_MAP.json` and `RELEASE_ASSETS.json`, and verify every selected input's bytes and SHA-256. The interrupted first build was not a complete release.
2. Restore a new work tree using only repository files and release assets, without access to the original archive. Rerun the finite check there. The successful test above deliberately used the original root and therefore does not certify the reduced package's closure.
3. Test the scalar/resource/endpoint entry points and the relocated analytic entry points in that restored tree. The mass and selected-moment adapters intentionally stop before obsolete endpoints; their output scopes must retain that distinction.
4. State Linux/WSL requirements for the parallel modular verifier (`fcntl`, `fork`, GMP/FLINT builds), Python version, `python-flint==0.9.0`, and that assertions must remain enabled. A syntax check is not an execution test.
5. Provide the exact commands for rebuilding the C++ programs and for fresh determinant recomputation. Saved-record checking must not be labelled determinant replay.
6. After upload, verify the public tagged commit and each release asset identity. Only then supply a public repository/release URL to the submission. A local directory or a partial ZIP is not that evidence.

The publisher should retain the unchanged manuscript hash receipt. This review is an independent agent audit of portability and evidence scope, not a formal proof-assistant verification or human mathematical certification.
