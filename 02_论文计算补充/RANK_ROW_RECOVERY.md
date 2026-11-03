# Recovery of complete source coefficients

The 32 entries in `indices/RECOVERED_RANK_ROWS.json` replace a request for
unavailable historical raw-JSONL bytes by an exact recovery of their mathematical
coefficient rows. Original matrix and certificate files remain byte-for-byte
unchanged. Historical raw filenames and SHA-256 values identify historical
files; they are not the names or hashes of the newly serialized rows.

Run from `02_论文计算补充` after obtaining its six original data volumes:

```
python tools/recover_rank_rows.py --output /absolute/path/to/new/recovered_rows
```

The tool reads authenticated matrix objects directly from `data/*.zip`, or from
the expanded `objects/` store. `--input-tree /absolute/path/to/input_tree` instead
reads the authenticated logical input paths. `--d 5` restricts the selection to
one indexed parameter. Every output directory must be new.

For a normalized row `sum_j c_j X^j`, `X=x^2`, the stored matrix consists of
`c_3,...,c_N`. The three identities proved in Theorem 6.1 give

```
c_0 = 0
c_2 = -sum(binomial(2*j,j)*c_j for j >= 3)/30
c_1 = 12*c_2.
```

Indeed, the `j=1,2` terms in the total-moment sum equal `2*12*c_2+6*c_2`.
The tool checks divisibility by 30 before the division. Multiplication by the
recorded `divided_small_content` recovers the full integer row. It also checks
that every prime divisor of that content is at most the strict characteristic
gate, that coefficients beyond the row's declared degree vanish, that every
row parameter `(d,t,m)` occurs exactly once in the specified family, and that
the parameters' monic filters agree with a fresh polynomial-lcm construction.
Certificate, matrix, and historical raw-reference identities are cross-bound.

These conditions recover and authenticate the candidate coefficients. The
low-moment identities do not prove that the row is the output of the original
Phi recurrence. The separate source check is the next step.

## Source recurrence interface

Add `--jobs` to export one `phi_d*.in` file per parameter:

```
python tools/recover_rank_rows.py --output /absolute/path/to/new/rows_and_jobs --jobs
```

Each stream contains **all** restored `(d,t,m)` rows and all their coefficients
from `X^0` to `X^t`, in the input format of the existing public source
`agents/mixed_source_independent/phi_raw_modular.cpp`
(SHA-256 `99083c4362ea0346e958b926cc731e30979fdcbe47e0e5a9cb073aa81c722bc6`).
After restoring the input tree, compile that source with GMP and run it on
each generated stream. For example:

```
g++ -O2 -std=c++17 /absolute/path/to/input_tree/agents/mixed_source_independent/phi_raw_modular.cpp -lgmpxx -lgmp -o /absolute/path/to/work/phi_raw_modular
/absolute/path/to/work/phi_raw_modular /absolute/path/to/new/rows_and_jobs/phi_d5.in
```

The recurrence uses `f_t=(-1)^t F_t`, which is the source-row sign convention
in the integer matrices. It independently constructs the source polynomial
by its quadratic power-sum recurrence and reduces it by the two monomial
rules for Phi. The filter input is rebuilt by this recovery tool from `d,t`.

To make a modular comparison imply equality over the integers, the job writer
uses a power of two strictly greater than twice the proved bound

```
H(t,g,m) = 9^(t+m)
           * sum((2*t/(t+i))*binomial(t+i,t-i)*13^i for i=0,...,t)
           * sum(abs(g_i)*9^i).
```

The power-sum formula gives the source coefficient majorant, while each Phi
rewrite reduces total degree and has absolute branching weight at most three.
Consequently the l1 norm of the reduction of a degree-D monomial is at most
`3^D`, giving weight nine in squared variables. The writer checks each term's
integrality and that every submitted coefficient satisfies this bound. Equal
residues then determine equal integer coefficients. Job generation records no
claim that the external recurrence has already been run.

The new recovery receipt explicitly records `historical_bytes_reconstructed:
false` and `source_recurrence_checked: false`. A complete source verification
requires successful external execution for every exported row, not merely
recovery success or a matching aggregate count.
