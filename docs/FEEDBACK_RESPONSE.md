# Supplementary response for the unchanged manuscript

This document accompanies *Explicit connectivity bounds for Markoff graphs*. The manuscript claims connectivity for every prime `5 <= p <= 10^66` and every prime `p >= 8*10^152`. Both ranges use computational certificates. The paper itself remains unchanged; its file identities are recorded separately in the repository.

## Computational supplement

The substantive issue raised by the feedback is that the PDF refers to source code and certificate files as part of its proof without supplying a retrieval location. A reproducible supplement must contain more than checking scripts: it must also contain the matrices and source rows, normalization/content data, selected minors and determinant certificates, exhaustive finite candidate list and compositeness witnesses, spectral node streams and interval certificates, and the exact endpoint inputs.

The release README and dependency index identify the fixed version, its data assets, file hashes, run commands, expected outputs and the mathematical statement checked by each command. The distinction between inspecting a saved certificate, checking its binding to a source, regenerating its rows, and recomputing its determinants is explicit. A source hash or archive-integrity check establishes identity, not the validity of a mathematical claim.

The main referee request also asks for a URL or DOI in Sections 1 and 10. Publishing the supplement and entering its fixed release URL in the submission's repository field supplies the missing external access. It does not insert those citations into the unchanged PDF. Any requirement to amend the PDF itself remains a task for a later manuscript version. The supplied decision text is insufficient to determine whether the site's versioning workflow accepts unchanged PDFs with new attachments or metadata.

## Earlier result of Martin [13]

Martin's [version 2, Theorem 1.1](https://arxiv.org/abs/2510.07577v2), proves the relevant orbit classification for primes `p` for which `24504480` does not divide `p^2-1`. In the classical case `kappa=0`, the exceptional orbit is the origin, so the theorem gives the nonorigin connectivity statement. The abstract and introduction state that the residue-class matrix criterion was checked for every prime power `d <= 17`.

The polynomial orbit-relation and residue-class rank method in the present manuscript follows Martin. The additional computation extends the certified parameter set and addresses the surviving finite candidates with specified conditional certificates. The frozen introduction cites the method but does not spell out Martin's earlier explicit modulus; this note records that quantitative comparison. The new numerical ranges should be assessed against Martin's result, rather than treating the rank method itself as new.

## Kingsbury-Neuschotz [10]: versions and title

The frozen bibliography explicitly cites [arXiv:2603.04096v1](https://arxiv.org/abs/2603.04096v1). Its official title is *Strong Approximation for the Character Variety of the Four-Times Punctured Sphere*, matching the version-specific title in the bibliography.

The current [version 3](https://arxiv.org/abs/2603.04096v3) is titled *Strong Approximation for the Relative Character Variety of the Four-Times Punctured Sphere*. Thus the feedback names the current title, while the bibliography names an earlier explicitly cited version. Version 3 also reports substantive proof corrections and a repaired gap. A future update should check the scope of the results cited before changing the version and title together. The reference supplies background here; the explicit prime ranges do not use its theorem as a proof premise.

## Numerical margins

[NUMERICAL_MARGINS.md](NUMERICAL_MARGINS.md) lists all twelve rows of Lemma 9.1 and the scalar gates (38), (40), and (41). The companion JSON records the source receipt hash and the exact scope of each calculation. The parameterized standard-library program `tools/numerical_margins.py` reproduces the table.

For the twelve endpoint rows, Decimal values only propose short rational brackets. The program then verifies those brackets with exact rational powers. The tightest stored upper expression is the count crossover `x_m`, leaving a rigorously certified margin of at least `0.907067` below `28885`.

The three original integer inequalities are also checked exactly. Their normalized margins are displayed separately: approximately `0.736546955` below the ordinary cap `8872`, `0.438484841` below the transport cutoff `14086`, and `52.792042050` between `28885` and the maximum constant allowed by the final scalar gate. Those three decimal displays are explanatory, not directed interval certificates. Changing an upstream constant requires rerunning the dependent inequalities; the displayed slack alone is not a sensitivity theorem.

## Computer assistance, process and verification notes

The abstract states the claimed ranges, while the introduction and keywords explicitly identify computer assistance. The submission's description or process field can also state that both ranges use computer-assisted proofs. Adding that phrase inside the abstract would require a later PDF revision.

The declared role `ai-discovered` should be accompanied by factual process and verification notes. These should identify the AI tools and research roles actually confirmed for this project, distinguish AI review from exact computational replay, list the computations genuinely rerun, and identify archived computations not independently regenerated in the current review. They should not infer a complete human verification or proof-assistant certification from an automated review, a successful build, file hashes, or a tagged release.

The feedback's observation that the theoretical reductions read showed no defect is a limited positive audit result. Verification of the missing calculations still requires the supplied inputs and the specified mathematical checking or regeneration procedures.
