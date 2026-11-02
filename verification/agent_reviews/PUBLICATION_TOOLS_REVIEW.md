# Independent review of publication tools and restored inputs

Date: 2026-11-02. This review concerns the computational supplement for the frozen manuscript. The reviewer did not modify the manuscript, create a repository, upload files, or certify the complete mathematical proof.

## Publication and restoration tools

The reviewed tools are `restore_inputs.py`, `verify_artifacts.py`, `run_checks.py`, `bundle_schema.py`, and the private publication launcher. The following safeguards were added and tested:

- Validate all logical paths, storage paths, hashes, byte counts, aliases, and release-asset declarations before restoration. Reject path traversal, drive-qualified or absolute paths, control characters, duplicate logical paths or asset names, alias length conflicts, undeclared assets, and unused asset declarations.
- Extract only the SHA-addressed ZIP members explicitly named in the input map. Hard links are confined to the new writable work tree; the repository and release payloads are not linked as writable work files.
- Use the same path validation in the artifact verifier. A work-only verification now checks restored release objects even when the ZIP cache is not supplied. Repository-only verification is labelled accordingly.
- Restrict credential-bearing publication requests and redirects to the HTTPS GitHub API/upload hosts. Credentials remain in memory and are not written to state files, URLs, or command arguments.
- Stage only files named and SHA-verified in `PUBLIC_FILES.json`, using literal Git path specifications. Unexpected staged files or an unexpected committed tree stop publication. Runtime caches and partial files are excluded from public-manifest generation.

Fourteen adversarial manifest/path/request tests passed. A separate work-only verifier regression test confirmed that a corrupted restored release object is rejected without an asset cache. The actual prepared input map (1,625 logical paths and 19 release assets at the time of this check) passed schema validation.

The prepared directory was scanned read-only for recognized GitHub/OpenAI access-token patterns, credential-bearing GitHub URLs, and private-key headers. No matches were found in 1,440 scanned non-cache files. Two generated Python cache files were identified as runtime artifacts and excluded from publication. This is a pattern-based check, not a universal guarantee about every possible confidential string.

## Tests using only the final repository and restored public inputs

These tests used the new `review_work` tree restored from the repository and its 19 release ZIPs. The original research archive root was not supplied to these commands.

| Stage | Actual result | Scope |
|---|---|---|
| `scalar` | PASS | Exact scalar/exponent comparisons in the supplied checker |
| `finite-records` | PASS | Fresh reconstruction of 27,985 CRT candidates, 26,210 recomputed compositeness witnesses, and source-bound saved records for 1,775 pointwise cases in 90 batches |
| `resource` | PASS | Exact integer/Fraction replay of 194 moment objects for 79 primes, including exponent tails and the resource bound |
| `endpoint` | PASS after the metadata adaptation below | Twelve endpoint rows and the final integer comparisons, using the unchanged integer/Fraction mathematical calculations |
| `small-graphs` | PASS | Fresh Windows C++ compilation and complete BFS for all 211 primes in the indexed range; output CSV compared line by line |

The finite-record test explicitly lists 92 omitted historical executable-byte identities: two source generators and 90 batch determinant executables. Their bound original C++ source hashes are checked. This test does **not** rebuild those executables or recompute the 1,775 large determinants or the 57 general rank certificates.

The initial endpoint test failed because the archived checker also hashed historical proof reports, old endpoint diagnostics, and handover metadata omitted from the reduced supplement. Inspection of its entire input-reading interface showed that none of these 21 missing references supplies numerical data to the current endpoint calculation. The portable adapter therefore records their exact historical path/hash declarations separately; it does not claim to verify those omitted bytes or the old audit/endpoint claims. It retains the three required control-file hashes/statuses, all remaining inherited-file identity checks, and every numerical check. The omission list is fixed and explicit in `indices/PORTABILITY_CHANGES.json` and in the fresh output.

All arithmetic helper functions are AST-identical to the archived endpoint checker. From `exact_budget_delta` through the twelve endpoint calculations and final integer comparison, the main-program source is byte-identical. The successful adapter reports 102 program gates and 12 endpoint rows; the gate count includes identity/declaration checks and should not be described as 102 distinct mathematical theorems.

## Evidence and remaining boundary

Local receipts are under `audit/fresh_public_replays/`; the repaired endpoint receipt is in `endpoint_repaired/endpoint/FRACTION_CONTINUOUS_BUDGET_REVIEW.json`. Supporting reviews are `PUBLICATION_TOOL_NEGATIVE_TESTS.json`, `WORK_ONLY_VERIFIER_NEGATIVE_TEST.json`, `ENDPOINT_PORTABILITY_OMISSIONS.json`, `ENDPOINT_MATH_BLOCK_UNCHANGED.json`, and `PUBLIC_STAGING_SCAN.json`.

These checks establish the stated portable execution and artifact-validation scopes. They do not replace the independent full determinant regenerations, every inherited transcendental-node computation, a human mathematical review, or formal proof-assistant verification. The public release must still be verified after upload against its fixed tag and all asset digests. A source-only Git clone does not contain the large release payloads.
