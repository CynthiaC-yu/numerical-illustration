# Numerical illustration of theoretical results

This repository provides numerical illustrations for *Feedback Stabilization for Sampled Linear Systems with Control-linear Noise*. It plots the observed Riccati convergence boundaries for five diagonal matrices `A_i` and two input cases, `B = I_5` and `B = ones(5,1)`. It also simulates the controlled system using the converged Riccati matrix.

## Run

The code has been run with MATLAB R2025a. Run the commands below from the project folder.

| Task | Entry point | Requirements |
| --- | --- | --- |
| Main boundary experiment and sample-path figures | `main.m` | MATLAB and Parallel Computing Toolbox |
| Four-tolerance comparison for all five matrices | `outputs/A1_A5_relTol_comparison/run_A1_A5_relTol_comparison.m` | MATLAB and Parallel Computing Toolbox |
| Redraw the saved comparison | `outputs/A1_A5_relTol_comparison/plot_A1_A5_relTol_comparison.m` | MATLAB |
| Reproduce and save the T=20 sample paths | `outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/run_sample_paths.m` | MATLAB |
| Optional high-precision boundary experiment | `main_vpa.m` | MATLAB, Parallel Computing Toolbox, and Symbolic Math Toolbox |

To run the main experiment:

```matlab
main
```

The script adds `dependencies/` to the path, starts or reuses a process pool, and computes 410 boundary searches: 255 for `B = I_5` and 155 for `B = ones(5,1)`. It saves both boundary plots and their data in a new output folder. Sections 9–10 then generate the two sample-path figures and leave their results in the MATLAB workspace.

To run only the sample paths, execute Sections 0–5, followed by Sections 9–10. The models, solver results, simulation settings, and trajectory data remain in `pathModels`, `pathProbes`, `pathCfg`, and `pathResults`.

A saved sample-path experiment with `relTol = 1e-6` is available in [`outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/`](outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/). To reproduce and save both sample-path figures directly:

```matlab
run(fullfile('outputs', 'sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6', 'run_sample_paths.m'));
```

The standalone script locates the project from its own folder and replaces the results and figures saved beside it.

## Boundary experiment

For `i = 1,...,5`,

```matlab
A_i = ((i-1)/5) * eye(5) + diag((0:4)/10);
```

The following settings are used by `main.m` and the saved main boundary experiment:

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
| Boundary search | Bisection in `sigma^2`, initially in `[0,4]` |
| Search precision | 2% relative bracket width in `sigma^2` |
| Continuous solver at `tau = 0` | `ode15s`, maximum time `5000`, `RelTol = 1e-9`, `AbsTol = 1e-11` |

For positive `tau`, the convergence criterion is

```matlab
norm(Wnew - Wold, 'fro') <= solverCfg.relTol * norm(Wold, 'fro')
```

The accepted `Wnew` is returned. At `tau = 0`, the solver instead checks a scaled Riccati derivative residual. Each boundary candidate starts from `P`.

The bisection stops when

```matlab
sigmaSquareHigh - sigmaSquareLow <= ...
    searchCfg.sigmaSquaredRelTol * max(sigmaSquareHigh, 1e-12)
```

The search retains a lower endpoint that passes the convergence check and an upper endpoint that fails it. It reports the lower variance endpoint as the observed `sigma_max^2`. The returned fields `sigmaSquareLow` and `sigmaSquareHigh` are both in variance units. Each midpoint solve receives `model.sigma = sqrt(sigmaSquareMid)`.

The scan stores `sqrt(sigmaSquareLow)` in `sigmaCurve`; the plotting function squares these values. A failed convergence check is a numerical outcome under the stated tolerances and limits. If the initial endpoints do not form a valid bracket, the search returns `NaN` endpoints.

The saved main boundary experiment uses the `epsilon = 1e-5` slice of the same completed rerun as the four-tolerance comparison below.

## Stopping-tolerance comparison for all five matrices

[`outputs/A1_A5_relTol_comparison/`](outputs/A1_A5_relTol_comparison/) contains both input cases for `A_1` through `A_5` at `epsilon = 1e-3, 1e-4, 1e-5, 1e-6`, where `epsilon` denotes `solverCfg.relTol`. All four settings use the relative update criterion above for positive `tau`. The model, sampling grids, quadrature, search precision, and solver limits are otherwise the same.

Each figure has five colors for the matrices and four line styles for the tolerances. The two groups have separate legend boxes. Diamonds, circles, squares, and triangles identify the four tolerances in the order listed above. Their markers are staggered along the common sampling grid to distinguish nearly overlapping curves.

The comparison run script reads the model and settings from its adjacent `results.mat`, reruns all 1,640 boundary searches, and replaces that file and both figures:

```matlab
addpath(fullfile('outputs', 'A1_A5_relTol_comparison'), '-begin');
comparison = run_A1_A5_relTol_comparison;
```

To redraw both figures from the saved comparison:

```matlab
run(fullfile('outputs', 'A1_A5_relTol_comparison', 'plot_A1_A5_relTol_comparison.m'));
```

The MAT file stores the settings, boundary diagnostics, and curve arrays in the order `(tau, A_i, epsilon)`. The dimensions are `51-by-5-by-4` for `cases(1)`, corresponding to `B = I_5`, and `31-by-5-by-4` for `cases(2)`, corresponding to `B = ones(5,1)`. The tolerance order is `[1e-3, 1e-4, 1e-5, 1e-6]`.

In each case, `sigmaSquared` stores the lower variance endpoints directly and `sigmaCurve` stores their square roots.

The estimates differ at 58 of 410 sampled `(A_i, B, tau)` combinations. The four estimates agree exactly at the remaining 352 combinations. The largest observed relative difference is 2.04%, using `epsilon = 1e-6` as the reference and unrounded values for the calculation.

| Input case | Sampled combinations | Combinations with different estimates | Maximum relative difference |
| --- | ---: | ---: | ---: |
| `B = I_5` | 255 | 34 | 2.04% |
| `B = ones(5,1)` | 155 | 24 | 1.89% |

For each sampled combination, let `v` be its four variance estimates. The reported difference is

```matlab
100 * max(abs(v - v(end))) / v(end)
```

The 2% search-width tolerance uses the upper bracket endpoint as its denominator. The comparison uses the `epsilon = 1e-6` estimate as its denominator.

Each case stores a `changedRows` table with one row per combination where the estimates differ, including the four estimates and the maximum relative difference as a percentage:

```matlab
comparison = load(fullfile('outputs', 'A1_A5_relTol_comparison', 'results.mat'));
comparison.cases(1).changedRows
comparison.cases(2).changedRows
```

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

## Optional high-precision workflow

The saved boundary and comparison experiments above use double-precision arithmetic. The additional `main_vpa.m` workflow uses 50-digit VPA arithmetic and the following convergence criterion:

```matlab
norm(Wnew - Wold, 'fro') <= 1e-3 * min(1, norm(Wold, 'fro'))
```

The model matrices, sampling grids, quadrature nodes, search precision, iteration limit, and matrix-norm limit match the boundary experiment above. At `tau = 0`, `integrate_riccati_vpa.m` uses adaptive Dormand–Prince 5(4) integration in VPA, with `odeMaxStep = 0.1`, `odeMaxTime = 5000`, `odeRelTol = 1e-9`, and `odeAbsTol = 1e-11`.

Set `nWorkers` in `main_vpa.m` to choose the process-pool size; `[]` uses MATLAB's default or an existing pool. This workflow scans `B = ones(5,1)` first and saves each completed boundary. Restarting the script starts a new run; checkpoints do not provide automatic resumption.

## Dependencies

| File | Role |
| --- | --- |
| `riccati.m` | Solves one `(A, B, tau, sigma)` from `model.P`; supports single and multiple inputs. Uses a fixed-point iteration for positive `tau` and a continuous Riccati ODE at zero. |
| `find_sigma_boundary.m` | Bisects `sigma^2` for one `(A, B, tau)` and retains both endpoints and solver results. Returns `NaN` estimates if the initial bracket is invalid. |
| `scan_sigma_curves.m` | Runs independent `(A_i, tau)` boundary searches with `parfor` and reports progress. |
| `plot_sigma_curves.m` | Squares the estimated `sigma` values and plots the boundary curves on a logarithmic vertical axis. |
| `sample_path.m` | Computes the feedback gain from a supplied `W`, simulates one selected model, and returns its figure and trajectory data. |
| `riccati_vpa.m` | High-precision Riccati solver with the min criterion. |
| `integrate_riccati_vpa.m` | High-precision continuous Riccati integration at `tau = 0`. |
| `find_sigma_boundary_vpa.m` | Bisects `sigma^2` using the VPA solver. |
| `scan_sigma_curves_vpa.m` | Parallel VPA boundary scans with progress reports and per-boundary checkpoints. |

## Saved outputs

| Folder under `outputs/` | Contents and role |
| --- | --- |
| [`boundary_sigmaSquared_relTol1e-5/`](outputs/boundary_sigmaSquared_relTol1e-5/) | Main boundary experiment at `relTol = 1e-5`: `results.mat`, `B_I5.fig/png`, and `B_ones5.fig/png`. |
| [`A1_A5_relTol_comparison/`](outputs/A1_A5_relTol_comparison/) | Four-tolerance comparison for all five matrices: run and plotting scripts, `results.mat`, and `B_I5_A1_A5_relTol_comparison.fig/png` and `B_ones5_A1_A5_relTol_comparison.fig/png`. |
| [`sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/`](outputs/sample_paths_A3_tau0p2_sigma0p1_T20_relTol1e-6/) | Additional sample-path example: `run_sample_paths.m`, `results.mat`, and `A3_B_I5.fig/png` and `A3_B_ones5.fig/png`. |

Both boundary experiments use direct bisection in `sigma^2`. Their MAT files record `searchVariable = 'sigmaSquared'`.

In the saved main boundary experiment, `sigmaCurveI` and `sigmaCurve` store sigma in `51-by-5` and `31-by-5` arrays, with one row per sampling interval and one column per matrix. Each MAT file also contains the model matrices, sampling grids, solver and search settings, convergence criterion, boundary diagnostics, and scan times.

The standalone sample-path script saves both models, all solver and simulation settings, the converged matrices, gains, and trajectory data. PNG files use 300 dpi.

New `main.m` runs save their boundary data and both FIG/PNG pairs under `outputs/<timestamp>_relative_tol_<relTol>/`.

`main_vpa.m` writes to `outputs/<timestamp>_vpa/`, saving startup settings, per-boundary checkpoints, case results, and both boundary figures. For runs with saved configurations, use those settings when comparing results.
