# Structure Transport in Factorized Controller Synthesis under Linear Equality Constraints

MATLAB code for the numerical examples in the manuscript:

"Structure Transport in Factorized Controller Synthesis under Linear Equality Constraints."

The scripts reproduce the numerical results, screening tests, and figures reported in the paper.

## Requirements

The code was tested with:

- MATLAB R2024a
- YALMIP 20200116
- SDPT3 4.0

Timing results reported in the paper were obtained on an Intel(R) Core(TM) Ultra 9 285H CPU at 2.90 GHz.

Before running the scripts, make sure YALMIP and SDPT3 are on the MATLAB path.

## Files

### `example1_identical_block.m`

Reproduces the identical-block H2 benchmark.

The script:

- evaluates the self-transport recovery-core certificate;
- verifies infeasibility of the self-transport slice for all `alpha > 0`;
- solves the non-self transport slice on the 101-point logarithmic grid `logspace(-2,2,101)`;
- reports the best tested certified bound, realized H2 cost, closed-loop spectral abscissa, and recovered core gain;
- recomputes the realized H2 cost of the external benchmark controller.

Representative values reported in the paper are

```text
alpha_best          = 0.9120
gamma_best          = 126.4345
realized H2^2       = 28.8542
spectral abscissa   = -0.4890
```

### `example1_timing.m`

Reproduces the computational comparison for the self-transport slice.

It compares:

- direct search: 101 complete H2 synthesis SDPs;
- proposed screening: one recovery-core certificate SDP.

On the machine listed above, the run reported in the paper gave approximately

```text
Direct grid search  : 11.4623 s
Proposed screening  : 0.1340 s
```

### `example2_mixed_affine.m`

Reproduces the mixed affine-equality example.

The controller satisfies

```text
k11 = 0
k23 = 1
k12 = k21
k22 + 2*k13 = 3
```

The script:

- constructs the prescribed transport-reference library;
- checks denominator compatibility;
- solves the recovery-core certificate SDP for every candidate;
- solves all surviving slices on the common 101-point logarithmic alpha grid;
- prints the best tested result for each reference family;
- reports the controller corresponding to the smallest certified bound;
- identifies the slice with the smallest realized H2 cost;

Representative results are

| Reference | tau | alpha  | gamma  | realized H2^2 |
| --------- | ---:| ------:| ------:| -------------:|
| `F1 = I3` | -   | 0.3020 | 2.2813 | 1.7314        |
| `F2 = Fb` | -   | 0.1738 | 2.0817 | 1.6747        |
| `F12+`    | 0.2 | 0.2089 | 2.1608 | 1.6684        |
| `F12-`    | 0.5 | 0.1738 | 2.0319 | 1.6422        |
| `F23+`    | 0.2 | 0.1585 | 2.1212 | 1.6851        |
| `F23-`    | 0.5 | 0.2291 | 2.0258 | 1.6562        |

## Notes on reproducibility

The numerical values can vary slightly with solver tolerances, operating system, and software versions. 

The files example1_corollary1.m and example2_corollary1.m provide independent implementations of the fixed-coefficient equalities of Corollary 1. The former uses nullspace elimination by default, while the latter imposes the equalities directly in YALMIP. Both yield results consistent with the explicit fixed-slice parameterizations.
