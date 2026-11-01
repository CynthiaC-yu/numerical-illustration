%% Numerical illustration of theoretical results
% Figures use the paper's visual style.

%% 0. Set up the experiment path
clear;
clc;

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'dependencies'), '-begin');

%% 1. Define the five A matrices
n = 5;
Aset = cell(5, 1);
for i = 1:5
    Aset{i} = ((i-1)/5) * eye(5) + diag(((0:(n-1))/10));
end

%% 2. Define input and cost matrices
Bset = cell(2, 1);
Bset{1} = eye(n);
Bset{2} = ones(n, 1);

% R is the identity in the corresponding input dimension.
Rset = cell(2, 1);
for i = 1:size(Bset, 1)
    m = size(Bset{i}, 2);
    Rset{i} = eye(m);
end

L = eye(n);

%% 3. Define tau grids and search precision
% Include tau=0; riccati uses its continuous ODE branch at that point.
tauSet = cell(2, 1);
tauSet{1} = (0:0.2:10).';
tauSet{2} = (0:0.2:6).';

% Bisect sigma; measure the relative bracket width in sigma^2.
searchCfg.sigmaSquaredRelTol = 0.02;
% Number of trapezoidal quadrature nodes for positive tau.
solverCfg.Ns = 201;

%% 4. Prepare the model for B = ones(5,1)
% The scan sets A and tau; the boundary search sets sigma for each candidate.
modelOnes.A = Aset{1};
modelOnes.B = Bset{2};
modelOnes.sigma = 0;
modelOnes.L = L;
modelOnes.R = Rset{2};
modelOnes.tau = tauSet{2}(2);
modelOnes.P = zeros(n); % P = 0

%% 5. Configure the Riccati solver
solverCfg.maxIter = 2e5;
solverCfg.maxWNorm = 1e12;
solverCfg.relTol = 1e-5;
solverCfg.progressInterval = 5; % Seconds between parallel progress reports.

% ODE integration settings for tau=0.
solverCfg.odeSolver = 'ode15s';
solverCfg.odeMaxTime = 5000;
solverCfg.odeRelTol = 1e-9;
solverCfg.odeAbsTol = 1e-11;

%% 6. Scan all A matrices for both input cases
% B = I_5
fprintf('\n=== B = I_5 ===\n');
modelI = modelOnes;
modelI.B = Bset{1};
modelI.R = Rset{1};

[sigmaCurveI, curveBoundaryI, curveSecondsI] = ...
    scan_sigma_curves(Aset, tauSet{1}, modelI, solverCfg, searchCfg);

% B = ones(5,1)
fprintf('\n=== B = ones(5,1) ===\n');
[sigmaCurve, curveBoundary, curveSeconds] = ...
    scan_sigma_curves(Aset, tauSet{2}, modelOnes, solverCfg, searchCfg);

%% 7. Plot both input cases
% Figure 1: B = I_5
[figI, axI] = plot_sigma_curves( ...
    tauSet{1}, sigmaCurveI, ...
    'Observed maximum converged \sigma^2 for ({A}_i, {I})', ...
    2);
ylim(axI, [1e-4, 2]);

% Figure 2: B = ones(5,1)
[figOnes, axOnes] = plot_sigma_curves( ...
    tauSet{2}, sigmaCurve, ...
    'Observed maximum converged \sigma^2 for ({A}_i, 1)', ...
    1);

%% 8. Save numerical results and figure files
% Use a separate output folder for each run.
runStamp = char(datetime('now', 'TimeZone', 'America/Chicago', ...
    'Format', 'yyyyMMdd_HHmmss_SSS'));
outputDir = fullfile(projectDir, 'outputs', ...
    sprintf('%s_relative_tol_%g', runStamp, solverCfg.relTol));
mkdir(outputDir);

convergenceCriterion = "norm(Wnew-Wold,'fro') <= relTol*norm(Wold,'fro')";
save(fullfile(outputDir, 'results.mat'), ...
    'Aset', 'Bset', 'Rset', 'L', 'tauSet', 'solverCfg', 'searchCfg', ...
    'convergenceCriterion', 'sigmaCurveI', 'curveBoundaryI', 'curveSecondsI', ...
    'sigmaCurve', 'curveBoundary', 'curveSeconds');

savefig(figI, fullfile(outputDir, 'B_I5.fig'));
exportgraphics(figI, fullfile(outputDir, 'B_I5.png'), ...
    'Resolution', 300);

savefig(figOnes, fullfile(outputDir, 'B_ones5.fig'));
exportgraphics(figOnes, fullfile(outputDir, 'B_ones5.png'), ...
    'Resolution', 300);

fprintf('\nResults and figures saved to:\n%s\n', outputDir);


%% 9. Configure the two sample-path cases
pathAIndex = 3;
pathTau = 0.2;
pathSigma = 0.1;

pathModels = cell(numel(Bset), 1);
for bIndex = 1:numel(Bset)
    pathModel = modelOnes;
    pathModel.A = Aset{pathAIndex};
    pathModel.B = Bset{bIndex};
    pathModel.R = Rset{bIndex};
    pathModel.tau = pathTau;
    pathModel.sigma = pathSigma;
    pathModels{bIndex} = pathModel;
end

% Both input cases use the same initial state and random seed.
pathCfg.x0 = ones(n, 1) / sqrt(n);
pathCfg.rngSeed = 123;
pathCfg.nPaths = 10000;
pathCfg.nShownPaths = 10;
pathCfg.nSampleIntervals = 100;
pathCfg.nSubsteps = 200;
% Final time = 100 * 0.2 = 20; Euler-Maruyama step = 0.2 / 200 = 0.001.


%% 10. Solve W and plot sample paths for each chosen model
pathCaseNames = {sprintf('A%d_B_I5', pathAIndex); ...
    sprintf('A%d_B_ones5', pathAIndex)};
pathProbes = cell(size(pathModels));
pathResults = cell(size(pathModels));
pathFigures = gobjects(numel(pathModels), 1);
pathAxes = gobjects(numel(pathModels), 1);

for bIndex = 1:numel(pathModels)
    pathModel = pathModels{bIndex};
    fprintf('\n=== %s | tau = %.3g | sigma = %.3g ===\n', ...
        pathCaseNames{bIndex}, pathModel.tau, pathModel.sigma);

    pathProbe = riccati(pathModel, solverCfg);
    pathProbes{bIndex} = pathProbe;
    if ~pathProbe.converged
        error('main:PathRiccatiNotConverged', ...
            '%s: Riccati solve returned %s.', pathCaseNames{bIndex}, pathProbe.status);
    end

    [pathFigures(bIndex), pathAxes(bIndex), pathResults{bIndex}] = ...
        sample_path(pathModel, pathProbe.Wfinal, pathCfg);
    set(pathFigures(bIndex), 'Name', pathCaseNames{bIndex}, 'NumberTitle', 'off');
end
