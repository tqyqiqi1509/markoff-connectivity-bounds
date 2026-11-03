# Review corrections for aiXiv:2611.00017

These corrections apply to the supplement at commit
`1f9a6da729d4d93163ba270667378173120d860c`. The prime ranges and numerical
constants are unchanged. The corrected paper is identified by the PDF SHA-256
in `THEOREM_INPUTS.json`. The manuscript is distributed separately.

| Item | Supported correction |
|---|---|
| R1 | Exact original middle-interval window engine is included in the selected inputs. |
| R2 | 47 additional logical paths provide the missing parameter aliases, plans, general-rank sources, kernels and checkpoints. They require 36 new byte objects; nine other recovered objects already occur in the original volumes. Generated d131/d137 matrices have explicit builders and expected hashes. |
| R3 | `indices/RECOVERED_RANK_ROWS.json` and `tools/recover_rank_rows.py` recover the complete coefficients of 32 normalized-matrix families. Historical raw JSONL hashes remain historical identities; generated files receive their own hashes. |
| R4 | The paper defines basis coefficients through the even part and takes even parts in the spectral identity. Proposition 9.3 specifies real `k>0`. |
| E1 | The kernel generator docstring now states 15,104 cells on `[0,32]`; its executable mathematics is unchanged. |
| Candidate join | `check_finite.py` binds the 57 parameters and requires the complete ordered classification sequence to equal the reconstructed candidates. |

## Data and execution

The six original data volumes are unchanged. `data/paper-data-07-review-fixes.zip`
contains only the 36 additional objects. `INPUTS.json` now selects 5,934 logical
paths from 4,980 unique byte objects. The first-layer research archives are
unchanged and are no longer needed to recover the dependencies listed in
`indices/REVIEW_DEPENDENCIES.json`.

Obtain the added volume from the v3.1.0 release. `download_assets.py` restores
the original 27 v3 assets and this additional volume using `RELEASE_ASSETS.json`.
When using the complete repository tree with all assets, run the upload checker
before expanding data. Then, from this directory, use new output locations:

```
python tools/prepare_data.py
python tools/materialize.py /absolute/path/to/new/input_tree
python tools/check_finite.py --project-root /absolute/path/to/new/input_tree --output /absolute/path/to/new/BASE.json
python tools/replay_review_fixes.py --project-root /absolute/path/to/new/input_tree --output /absolute/path/to/new/MIDDLE --build-131-137
python tools/recover_rank_rows.py --input-tree /absolute/path/to/new/input_tree --output /absolute/path/to/new/RECOVERED_ROWS
python tools/finite_latest.py --project-root /absolute/path/to/new/input_tree --output /absolute/path/to/new/LATEST.json --reconstruct
```

The review replay checks all 89 supplied-file bindings and seven
parameter/source families, then 397 direct and 57 window branches and their
prime DAGs. The build option reconstructs d131/d137 matrices. Optional
`--case9002` requires Linux and python-flint 0.9.0, with output inside the
restored input tree. It checks one mixed certificate, including its integer
minors and gcd support. See `DEPENDENCY_RECOVERY.md` for the recovery roles.

The raw recovery was exercised on all 32 families (10,760 rows and 2,787,140
complete coefficients). It verifies filters, row keys, content, integrality
and coefficient reconstruction. Source legality is a separate check: use
`--jobs` and execute the public Phi recurrence as described in
`RANK_ROW_RECOVERY.md`. Recovery success alone is not a source-recurrence
certificate. The full original recurrence for all 32 families and all 59/176
rank determinant replays was not rerun in this correction round.

Existing analytic certificates are preserved. Their producer hashes identify
the original execution sources; the kernel generator's sole change is its
docstring. Current maintained-source hashes are in `TOOL_PROVENANCE.json`.
