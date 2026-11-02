# Numerical margins for the frozen paper

The twelve rows below bound the expression in Lemma 9.1 using the stored rational upper bounds for U, W and the trace. Each interval was checked afresh by exact rational exponentiation; it is not obtained by rounding a displayed decimal.

| Trace branch | Endpoint | Certified positive margin below 28885 |
|---|---|---|
| trivial | zero_limit | [1488331153/200000, 3720827883/500000] |
| trivial | cut_128 | [6929432303/1000000, 433089519/62500] |
| continuous_profile | cut_128 | [349041147/50000, 6980822941/1000000] |
| continuous_profile | horizontal_to_hyperbola | [6674325709/1000000, 667432571/100000] |
| continuous_profile | low_cap | [3116012771/500000, 6232025543/1000000] |
| continuous_profile | mass_crossover | [800699701/500000, 1601399403/1000000] |
| continuous_profile | cut_10000 | [683378603/1000000, 170844651/250000] |
| middle_41_100 | cut_10000 | [426182687/200000, 532728359/250000] |
| middle_41_100 | cut_14086 | [287656329/250000, 1150625317/1000000] |
| high_355_1000 | cut_14086 | [839509691/500000, 1679019383/1000000] |
| high_355_1000 | count_crossover | [907067/1000000, 226767/250000] |
| high_355_1000 | right_edge | [8581021/1000000, 4290511/500000] |

The tightest row is the count crossover x_m, with positive margin at least 0.907067.

## Scalar gates

All three original integer inequalities were checked afresh with integer exponentiation. The following decimals are explanatory displays, not directed interval certificates.

| Equation | Normalized interpretation | Displayed margin |
|---|---|---|
| (38) | 8872 minus the ordinary-cap scalar | 0.73654695535716903085342984347750046036115384961262551072661390008110843 |
| (40) | 14086 minus 14000(5/4)^(12533/459000) | 0.4384848413206312546982056503400853120069699881022196671819686670243480 |
| (41) | Largest admissible endpoint constant minus 28885 | 52.7920420503972847377741253676202733897286862484287624580274821602344634 |

These margins quantify the current bound. They do not by themselves permit individual certified input constants to change: sensitivity requires rerunning the dependent comparisons.

## Reproduce

From the repository root, after restoring the proof data to `work`:

```sh
python tools/numerical_margins.py --endpoint-receipt work/control/continuous_budget_fraction_root_replay/FRACTION_CONTINUOUS_BUDGET_REVIEW.json --output-dir results/margins
```

Use a new output directory and keep Python assertions enabled. This command uses only the standard library. The JSON identifies the exact input and script hashes. The twelve intervals describe the stored upper expressions, rather than two-sided bounds on the exact mathematical optimum.
