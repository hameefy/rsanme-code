# Experiments for the revised RSANME manuscript (rsanme_v2.tex)

All code is in MATLAB (developed for R2026a). Nothing in this folder has produced results for the paper yet.
Only a tiny smoke test (dimension 200-500, 2 repetitions) was run to check that the code executes; its
output was deleted. **You run the experiments.**

## Layout
- `solvers/`  rsanme.m, rsanmeq.m (RSANME-Q), dfsanme.m, pandsr.m, liu2022.m, modified_algorithm21.m,
  ISCGPM_solver.m (logistic regression only, unchanged from your earlier code), helpers.
  Common interface `[x, info] = solver(F, x0, opts)`, with budgets tol=1e-6, 1000 iterations, 2000 F-evaluations.
  Every evaluation of F is counted, including projection-step evaluations.
- `problems/` test problems P1-P5, starting points of the manuscript, training problems TR1-TR5 (tuning only).
- `drivers/`  the scripts below.

## Run order (from `code/`: `addpath('drivers')`)
1. `run_tuning`          tunes RSANME and RSANME-Q on TR1-TR5 (random search 300 configs + refinement,
                         validated on a disjoint suite) and writes `tuned_params.mat/.txt`. Run first.
2. `verify_strong_monotonicity`  checks the hypothesis (mu>0, a < mu/(L+r)) at the computed solutions of P1-P5.
3. `run_monotone`        main benchmark (180 instances x 6 solvers), CSV + performance profiles in `results/`.
4. `run_sensitivity_a`   sensitivity of RSANME to a and of RSANME-Q to sigma.
5. `run_logistic`        12 LIBSVM data sets, 12 repetitions, 6 solvers.
6. `run_signal`          sparse signal recovery: run, no positive result, REMOVED from the paper (kept only for the record).
Each driver accepts `true` for a quick test, e.g. `run_monotone(true)`.
`run_monotone` at n=150000 is the long job. Reported timings need an otherwise idle machine.

## Settings that need your confirmation
- Signal experiment: the manuscript does not give n, k, the number of non-zeros or the number of instances.
  Defaults in `run_signal.m`: n=4096, k=1024, 64 non-zeros, 10 instances, noise variance 1e-12, gamma=0.5.
- Competitors use their published parameters (listed in Section 4.1 of the paper). PAND-SR with its published
  eta_k = 0.99^k(100+||F0||^2) can blow up on the signal problem (the smoke test showed this); it is run as published.
  If you prefer, compete tuned versions as well (change `par` in the drivers; the paper text says only
  RSANME/RSANME-Q are tuned).
- Liu 2022 is run in the monotone benchmark only; ISCGPM only in logistic regression (as in the first version).
- Data folder for logistic regression: `Experiments/Logistic_Regression` (edit `dataDir` in `run_logistic.m`).

## After running
Send back `results/*.csv`, the figures and `tuned_params.txt`. The `\todo{...}` markers in `rsanme_v2.tex`
and the `[PENDING]` items in `response_to_reviewers.tex` are then filled from the numbers. Nothing numerical
is claimed in the manuscript until then.
