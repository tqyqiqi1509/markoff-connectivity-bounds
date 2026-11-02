# Ordinary moment constants 21718 and 135

The full normalized recursion and its integer replay are
`prior_restored/baragar_current_run/weighted_window_certificate/weighted_159.json`
and `replay_weighted_159.json`. The high-side ordinary constant is extracted from
support mask 63 (`C_all_two_upper`), not from the joint constant 22098 of that old
endpoint. The low support tail at q=31 is below 134.749609, hence below 135.
The 64 support recursions, exact contents and masks remain in the supplied files.

The proof uses the ordinary count constants from this record, not its obsolete
p>=10^159 endpoint. Its original integer replay recomputes both. The public
instructions must state that distinction when reporting this replay.

After restoring a writable working root, the exact standard-library replay is:

```sh
python prior_restored/baragar_current_run/weighted_window_certificate/replay_weighted_integer.py \
  --integrals prior_restored/baragar_current_run/audit/selberg_mellin/mellin_density16384.json \
  --recursion prior_restored/baragar_current_run/weighted_window_certificate/weighted_159.json \
  --prefixes prior_restored/baragar_current_run/weighted_window_certificate/prefixes.json \
  --output NEW_ORDINARY_REPLAY.json
```

This command uses all saved node bytes, verifies the streams and block/prefix
sums, recomputes roots, tails, harmonic cuts and every recurrence. It does not
independently regenerate every transcendental enclosure. The expected old replay
status is `PASSED`, with `integer_fraction_only: true`, and the recorded joint
constant 22098. The claimed ordinary constants are the high support row and the
q=31 tail, rather than that joint constant.

The input streams are all in the release/path map:

- `audit/selberg_mellin/mellin_density16384.nodes.u64.gz` under the restored
  `prior_restored/baragar_current_run/` root;
- `audit/interval_antichain/arb_178_full32768.nodes.u64.gz` under that same root;
- `audit/disjoint_support/support_169_32768.subset_nodes.u64.gz` under that root;
- the metadata files with matching stems, and weighted `prefixes.prefixes.json.gz`.

The high-level sources `certify_mellin_integrals.py`,
`independent_interval_arb.py`, `certify_disjoint_support.py`,
`build_weighted_prefixes.py` and `certify_weighted_recursion.py` are supplied.
For a fresh Mellin computation in a writable copy, the local source's interface is
`--density 16384 --lambdas 1,2,3,4,6,8,12,18 --output NEW.json`.
For fresh harmonic prefixes:

```sh
python prior_restored/baragar_current_run/weighted_window_certificate/build_weighted_prefixes.py \
  --input prior_restored/baragar_current_run/audit/selberg_mellin/mellin_density16384.json \
  --output NEW_PREFIXES.json
```

These source-level commands have been read from their actual argument parsers;
they are not claims that every expensive original computation was rerun during
the publication audit. The frozen source records retain old absolute execution
paths as provenance; explicit command-line inputs select the restored files.
