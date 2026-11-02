# rsanme-code

MATLAB code, tuned parameters and result files for the paper

> H. Mohammad, H. Choi, A. B. Abubakar, M. Abdullahi, M. S. Sarkinbai,
> *On a positive spectral residual algorithm with a norm descent line search for nonlinear monotone equations* (under revision).

The code solves monotone systems `F(x) = 0` with `RSANME` (norm descent line search) and `RSANME-Q`
(quadratic-decrease line search) and compares them with Modified Algorithm 2.1, Algorithm 2.1 of Liu et al.,
DFSANME, PAND-SR and (logistic regression only) ISCGPM.

## Layout
- `solvers/`   all solvers, common interface `[x, info] = solver(F, x0, opts)`
- `problems/`  test problems P1-P5, six starting points, and the disjoint training suite TR1-TR5
- `drivers/`   experiment scripts
- `results/`   CSV files and figures used in the paper
- `tuned_params.txt`, `tuned_params.mat`   parameters selected by `run_tuning`

## Reproducing the experiments
Run from MATLAB in `drivers/`, in this order:
`run_tuning`, `verify_strong_monotonicity`, `run_monotone`, `run_sensitivity_a`, `run_logistic`.
Each driver accepts `quick = true` for a small smoke test.
The performance profiles of the paper are drawn from `results/results_monotone.csv` by `python drivers/make_profiles.py`
(needs numpy, pandas, matplotlib); `run_monotone` also writes plainer versions through `perfprof.m`.
All runs in the paper used MATLAB R2026a (Update 5) on Windows 11, Intel Core Ultra 9 275HX, 24 cores, 63 GB RAM.

## Data
The logistic regression experiment uses the public LIBSVM data sets
(https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/). They are not included. Download them and set `dataDir`
in `drivers/run_logistic.m`.

## Note on the sparse-signal experiment
`drivers/run_signal.m` is included for completeness. No method gave a useful reconstruction in our runs, so the
experiment is not part of the paper and its results are not included.

## Licence
MIT, see `LICENSE`.
