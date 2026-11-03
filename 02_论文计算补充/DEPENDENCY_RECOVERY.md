# R1/R2 dependency repairs

R1 was confirmed by checking the second-layer input list: the required window engine was absent. Its original bytes were recovered from the fourth-round archive and match e808f193f6f0f1c4ac9d22d0226a31f1e988123154e3e5318263716805fc4fb1. The original wrapper was rerun on all three middle intervals: 397 direct rows and 57 window rows passed, with every input Pocklington-DAG node checked.

R2 was confirmed for the seven parameter aliases, three legacy plans, four d131/d137 raw sources, d167/d173 extension rows, factored kernels and block certificates/checkpoints, the final d121 CRT checkpoint and the separate d163 parallel verifier. indices/REVIEW_RECOVERY_PROVENANCE.json gives exact archive/member identities for every added logical input. Historical prose is marked historical-only in indices/REVIEW_DEPENDENCIES.json; it is retained to satisfy frozen legacy execution hash gates, not as a mathematical premise.

The old build_full_matrix.py was actually executed for d131 and d137. Its output hashes equal the existing certified hashes. These large intermediate matrices are generated, not duplicated in the added data volume. The 7 parameter sidecars were checked against all source-row d,n,m keys and counts. Case 9002 was completely rerun under an isolated Linux Python 3.14.4 / python-flint 0.9.0 environment: both top inverse identities, compression, two maximal-minor modular reconstructions and the final smooth gcd check passed. The final COMPRESSED_MIXED_REVIEW.json exists. This is one mixed case, not all 176.

The maintained replay_review_fixes.py defaults to exact dependency binding plus the three middle arithmetic replays. --build-131-137 regenerates the two specified intermediate matrices. --case9002 requires Linux/python-flint and an output directory inside the restored project tree, because the frozen checker records relative paths.

No original archive, certificate or frozen source was edited. No claim is made that all 59 general and 176 mixed rank determinants were freshly regenerated in this repair. Restoring exact raw coefficients and checking identities is distinct from source legitimacy under the Phi recurrence.
