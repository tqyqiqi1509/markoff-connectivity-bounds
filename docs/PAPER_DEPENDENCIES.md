# Frozen paper to computational inputs

The paper-file SHA-256 values are in `indices/FROZEN_PAPER.json`. Original
logical paths below are resolved by `data/PATH_MAP.json`, including aliases.
The original JSON/certificate bytes are not edited to replace historical paths.

| Paper input | Entry points and exact input identities |
|---|---|
| Theorems 5.2/6.1, 57 general parameters | `indices/general_rank_inputs.json`, all 57 certificate paths, source/matrix rows and chosen-minor/journal dependencies |
| Theorem 1.1, finite residue classes | `tools/check_finite.py`; `agents/finite_extension_round20/CRT57_1E66_INDEPENDENT_CLASSIFICATION.json`; 57 congruence moduli |
| 26,210 compositeness witnesses | `records` in the final classification; exact proper factors or failing Miller-Rabin power chains |
| 1,775 conditional minors | `indices/POINTWISE_BATCHES.json`; 90 batch manifests, source rows for d=167/173, GMP and FLINT records and source programs |
| Small primes | `tools/small_graphs.cpp`, `indices/small_graphs.csv`; exact BFS through 1304 |
| Nine support parameters behind (18) | General certificates d=22,33,15,35,39,34,38,46,58; `agents/support_coupling_round7/certify_support_coupling.py` and its inputs |
| Proposition 7.1: 21718,135 | Ordinary spectral streams and weighted recursion under `prior_restored/baragar_current_run/`; see `ORDINARY_MOMENTS.md` |
| Proposition 7.1: 34195,28924 | `prior_restored/latest/positive157/complete/`, positive integer/node replays, `agents/analytic_mass_round3/MASS_MOMENT_CERTIFICATE.json` |
| Proposition 7.1: 159842,662,26/25 and tau constant | Exact Euler factors in `agents/joint_analytic/HALF_MODULUS_ENDPOINT_CERTIFICATE.json`, support-coupling inputs and endpoint/scalar checks |
| Lemma 8.1: 43540 | `agents/probability_mechanism_round6/SELECTED_MOMENT_CERTIFICATE.json`; `tools/verify_selected_moment_independent.py` |
| Lemma 8.3: resource constant | `tools/resource/RESOURCE_CERTIFICATE.json`, `verify_resource_rational.py`; 194 finite resource objects |
| Proposition 8.4: 177/500 and 2/5 | `agents/gaussian_fixed_tail_certificate_round34/EULER_FREQUENCY_CERTIFICATE.json`, resource certificates, portable frequency/variance checkers |
| Lemma 9.1: 12 endpoints; (38),(40),(41) | `tools/replay_continuous_budget_fraction.py`; `indices/ENDPOINT_RECEIPT.json`; `tools/numerical_margins.py` |

Every included input has a file digest and a selection reason. Excluded
historic narrative, exploratory progress, compiled executables, duplicate
payload bytes and third-party literature are not additional theorem premises.
Some immutable records mention these historical objects; those references are
provenance, not promises that a discarded research tree is reproduced.
