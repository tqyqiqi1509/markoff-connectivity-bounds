# Markoff connectivity: proof materials v3.1.1

Computational support for connectivity of the nonorigin solutions of
x² + y² + z² = xyz over F_p under exactly the three Vieta involutions,
for every prime **5 ≤ p ≤ 10^84** or **p ≥ 10^151**.

This repository contains two layers. The manuscript itself is excluded.

| Layer | Contents | Entry |
|---|---|---|
| 1 | Five archived research rounds, 21 original ZIP volumes | [Archive index](01_五轮总材料/ARCHIVES.json) |
| 2 | Current computation tools, exact inputs and certificates; six original data ZIP volumes plus one correction volume | [Proof-support guide](02_论文计算补充/README.md) |

**[Download the current v3.1 release](https://github.com/tqyqiqi1509/markoff-connectivity-bounds/releases/tag/v3.1.0).**
Large ZIP files are release attachments. Their original relative paths,
sizes and SHA-256 values are in [RELEASE_ASSETS.json](RELEASE_ASSETS.json).
The complete upload-file inventory is [UPLOAD_FILES.json](UPLOAD_FILES.json).

## Download and verify

Clone or download this repository, then run:

```
python download_assets.py
python 02_论文计算补充/tools/check_upload.py
```

The downloader restores all listed ZIP attachments to their intended paths and
checks their SHA-256 values. Then follow [the reproduction guide](02_论文计算补充/REPRODUCE.md); expand the proof-data volumes before restoring the readable input tree.
To download only the current paper support, run `python download_assets.py --layer 2`; full-bundle checking requires both layers.

## Current evidence

- 59 general-parameter rank certificates; 176 mixed rules, including all 165 triples from the eleven-element pool.
- 90 batches for 1775 base-range candidates, with the previously omitted external data restored.
- Fresh complete CRT/witness checks for the base and last finite intervals.
- Fresh generation of all 15104 frequency cells on [0,32], signed-trace checks, and the twelve final endpoint comparisons.

The exact theorem-to-input map is [THEOREM_INPUTS.json](02_论文计算补充/THEOREM_INPUTS.json).
The maintained tools and receipts distinguish fresh numerical checks from inherited deterministic rank certificates. The full Baragar conjecture still has the middle prime interval unresolved.

This v3.1 update preserves the original v3 release assets and Git history.

## Review corrections

This update adds the corrections described in [REVIEW_FIXES.md](02_论文计算补充/REVIEW_FIXES.md). The 27 original v3 assets are unchanged. The additional `paper-data-07-review-fixes.zip` is attached to v3.1.0. The downloader obtains all 28 required ZIPs across the two releases. The prime ranges and numerical constants are unchanged.
