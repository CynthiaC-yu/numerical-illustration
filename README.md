# Numerical illustration of theoretical results

This repository provides numerical illustrations for *Feedback Stabilization for Sampled Linear Systems with Control-linear Noise*. It plots the observed Riccati convergence boundaries for five diagonal matrices `A_i` and two input cases, `B = I_5` and `B = ones(5,1)`. It also simulates the controlled system using the converged Riccati matrix.

## Run

The main workflow uses MATLAB and Parallel Computing Toolbox. It has been run with MATLAB R2025a. Open the project folder in MATLAB and run:

```matlab
main
```

The script adds `dependencies/` to the path, starts or reuses a process pool, computes both boundary plots, and saves their data and figures. Sections 9–10 then generate the two sample-path figures.

To run only the sample paths, execute Sections 0–5, followed by Sections 9–10. The models, solver results, simulation settings, and trajectory data remain in `pathModels`, `pathProbes`, `pathCfg`, and `pathResults`.

A saved sample-path experiment with `relTol = 1e-6` is available in [`outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/`](outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/). Run its `run_sample_paths.m` to reproduce and save both sample-path figures without running the boundary scans. The script locates the project from its own folder.

## Boundary experiment

For `i = 1,...,5`,

```matlab
A_i = ((i-1)/5) * eye(5) + diag((0:4)/10);
```

The following settings are used by the current `main.m`:

| Setting | Value |
| --- | --- |
| Input matrices | `B = eye(5)` and `B = ones(5,1)` |
| Cost matrices | `L = eye(5)`, `R = eye(size(B,2))` |
| Initial Riccati matrix | `P = zeros(5)` |
| Sampling intervals | `0:0.2:10` for `B = I_5`; `0:0.2:6` for `B = ones(5,1)` |
| Matrix integral | 201-point trapezoidal quadrature |
| Fixed-point tolerance | `relTol = 1e-5` |
| Maximum sampled-map updates | `maxIter = 2e5` |
| Frobenius norm limit | `maxWNorm = 1e12` |
| Boundary search | Bisection in `sigma`, initially in `[0,2]` |
| Search precision | 2% relative bracket width in `sigma^2` |
| Continuous solver at `tau = 0` | `ode15s`, maximum time `5000`, `RelTol = 1e-9`, `AbsTol = 1e-11` |

For positive `tau`, the convergence criterion is

```matlab
norm(Wnew - Wold, 'fro') <= solverCfg.relTol * norm(Wold, 'fro')
```

The accepted `Wnew` is returned. At `tau = 0`, the solver instead checks a scaled Riccati derivative residual. Each boundary candidate starts from `P`.

The bisection stops when

```matlab
sigmaHigh^2 - sigmaLow^2 <= ...
    searchCfg.sigmaSquaredRelTol * max(sigmaHigh^2, 1e-12)
```

The curve uses `sigmaLow`, the tested convergent endpoint. A failed convergence check is a numerical outcome under the stated tolerances and limits.

## Stopping-tolerance comparison for A_3

[`outputs/A3_relTol_comparison/`](outputs/A3_relTol_comparison/) contains the comparison for `A = A_3` and both input cases. The legend uses `epsilon` for `solverCfg.relTol`, with values `1e-3`, `1e-4`, `1e-5`, and `1e-6`. For positive `tau`, all four runs use the relative update criterion above. The model, sampling grids, quadrature, search precision, and solver limits are otherwise the same.

To rerun all 328 boundary searches and save the results and both figures:

```matlab
addpath(fullfile('outputs', 'A3_relTol_comparison'), '-begin');
comparison = run_A3_relTol_comparison;
```

To redraw the figures from the saved data:

```matlab
run(fullfile('outputs', 'A3_relTol_comparison', 'plot_A3_relTol_comparison.m'));
```

Both entries use the adjacent `results.mat`. Diamonds, circles, squares, and triangles identify the four tolerances in the order listed above. Within each figure, all four curves use the same `tau` grid; their markers are staggered to distinguish nearly overlapping lines.

The table lists all 7 of the 82 sampled `(B, tau)` pairs where the results differ. Values of `sigma_max^2` are displayed to three significant figures. Maximum relative differences are computed from unrounded values against the `epsilon = 1e-6` result. All other sampled pairs give identical values.

| B | tau | epsilon=1e-3 | epsilon=1e-4 | epsilon=1e-5 | epsilon=1e-6 | Max. relative difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `I_5` | 0.400 | 0.446 | 0.446 | 0.446 | 0.441 | 1.18% |
| `I_5` | 2.60 | 0.0413 | 0.0405 | 0.0405 | 0.0405 | 1.95% |
| `I_5` | 4.80 | 0.00222 | 0.00220 | 0.00220 | 0.00220 | 1.04% |
| `I_5` | 8.00 | 2.21e-5 | 2.21e-5 | 2.18e-5 | 2.18e-5 | 1.31% |
| `ones(5,1)` | 3.80 | 2.66e-5 | 2.63e-5 | 2.63e-5 | 2.63e-5 | 1.19% |
| `ones(5,1)` | 4.60 | 7.97e-6 | 7.97e-6 | 7.88e-6 | 7.88e-6 | 1.09% |
| `ones(5,1)` | 5.80 | 8.66e-7 | 8.52e-7 | 8.52e-7 | 8.52e-7 | 1.66% |

## Sample paths

Both cases use `A = A_3`, `L = eye(5)`, `P = zeros(5)`, `tau = 0.2`, and `sigma = 0.1`. Their input and control-cost matrices are `B = eye(5), R = eye(5)` and `B = ones(5,1), R = 1`.

| Setting | Value |
| --- | --- |
| Initial state | `x0 = ones(5,1)/sqrt(5)` |
| Riccati solver | `riccati.m`; `relTol = 1e-5` in `main.m`, `1e-6` in the saved example |
| Time integration | Euler–Maruyama, `dt = 0.001` |
| Simulation length | 100 sampling intervals; `T = 20` |
| Monte Carlo paths | 10000; 10 individual paths shown |
| Random seed | `123`, using `twister` in both cases |

The feedback is `u(k*tau+s) = K(s)*x(k*tau)`. Each path uses one scalar Wiener process. The figures show individual squared state norms, their Monte Carlo mean, and the mean at sampling instants on a logarithmic vertical axis.

An additional [`B = ones(5,1)`, `T = 100` figure](outputs/sample_path_A3_B_ones5_tau0p2_sigma0p1_T100_relTol1e-6/) is saved in PNG and FIG formats. It reuses the saved `relTol = 1e-6` matrix and the same initial state, random seed, path count, and step size, with `nSampleIntervals = 500`.

## Optional high-precision workflow

`main_vpa.m` additionally requires Symbolic Math Toolbox. It uses 50-digit VPA arithmetic and a different convergence criterion:

```matlab
norm(Wnew - Wold, 'fro') <= 1e-3 * min(1, norm(Wold, 'fro'))
```

The model matrices, sampling grids, quadrature nodes, search precision, iteration limit, and matrix-norm limit match the boundary experiment above. At `tau = 0`, `integrate_riccati_vpa.m` uses adaptive Dormand–Prince 5(4) integration in VPA, with `odeMaxStep = 0.1`, `odeMaxTime = 5000`, `odeRelTol = 1e-9`, and `odeAbsTol = 1e-11`.

Set `nWorkers` in `main_vpa.m` to choose the process-pool size; `[]` uses MATLAB's default or an existing pool. This workflow scans `B = ones(5,1)` first and saves each completed boundary. Restarting the script starts a new run; checkpoints do not provide automatic resumption.

## Dependencies

| File | Role |
| --- | --- |
| `riccati.m` | Solves one `(A, B, tau, sigma)` from `model.P`; supports single and multiple inputs. Uses a fixed-point iteration for positive `tau` and a continuous Riccati ODE at zero. |
| `find_sigma_boundary.m` | Bisects `sigma` for one `(A, B, tau)` and retains both endpoints and solver results. Returns `NaN` estimates if the initial bracket is invalid. |
| `scan_sigma_curves.m` | Runs independent `(A_i, tau)` boundary searches with `parfor` and reports progress. |
| `plot_sigma_curves.m` | Squares the estimated `sigma` values and plots the boundary curves on a logarithmic vertical axis. |
| `sample_path.m` | Computes the feedback gain from a supplied `W`, simulates one selected model, and returns its figure and trajectory data. |
| `riccati_vpa.m` | High-precision Riccati solver with the min criterion. |
| `integrate_riccati_vpa.m` | High-precision continuous Riccati integration at `tau = 0`. |
| `find_sigma_boundary_vpa.m` | Boundary search using the VPA solver. |
| `scan_sigma_curves_vpa.m` | Parallel VPA boundary scans with progress reports and per-boundary checkpoints. |

## Saved outputs

`main.m` saves the boundary experiment under `outputs/<timestamp>_relative_tol_<relTol>/`: `results.mat` and the `B_I5` and `B_ones5` figures in FIG and PNG formats. The MAT file contains the model matrices, sampling grids, solver and search settings, convergence criterion, curves, boundary diagnostics, and scan times.

`sigmaCurveI` and `sigmaCurve` store **sigma**, with one row per `tau` and one column per `A_i`. The plotting function performs the squaring.

The standalone sample-path script saves `results.mat`, `A3_B_I5.fig/png`, and `A3_B_ones5.fig/png` beside itself. Its MAT file includes both models, all solver and simulation settings, the converged matrices, gains, and stored trajectory data. PNG files use 300 dpi and compact outer margins.

`main_vpa.m` writes to `outputs/<timestamp>_vpa/`, saving startup settings, per-boundary checkpoints, case results, and both boundary figures. For runs with saved configurations, use those settings when comparing results.
