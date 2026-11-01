%% Sample paths for A_3, tau=0.2, sigma=0.1, and T=20
% Run this file to reproduce both cases using the project dependencies.
% Results and figures are saved beside this script.
clear;
clc;

outputDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(fileparts(outputDir));
addpath(fullfile(projectDir, 'dependencies'), '-begin');

%% Model and solver configuration
n = 5;
A = 0.4 * eye(n) + diag((0:(n-1)) / 10);
Bset = {eye(n); ones(n, 1)};
pathCaseNames = {'A3_B_I5'; 'A3_B_ones5'};
pathModels = cell(2, 1);
for bIndex = 1:2
    pathModels{bIndex} = struct('A', A, 'B', Bset{bIndex}, ...
        'L', eye(n), 'R', eye(size(Bset{bIndex}, 2)), ...
        'P', zeros(n), 'tau', 0.2, 'sigma', 0.1);
end

solverCfg.Ns = 201;
solverCfg.maxIter = 2e5;
solverCfg.maxWNorm = 1e12;
solverCfg.relTol = 1e-4;
solverCfg.progressInterval = 5;
% Retained from main.m; inactive for these positive-tau cases.
solverCfg.odeSolver = 'ode15s';
solverCfg.odeMaxTime = 5000;
solverCfg.odeRelTol = 1e-9;
solverCfg.odeAbsTol = 1e-11;

%% Simulation configuration
pathCfg.x0 = ones(n, 1) / sqrt(n);
pathCfg.rngSeed = 123;
pathCfg.nPaths = 10000;
pathCfg.nShownPaths = 10;
pathCfg.nSampleIntervals = 100;
pathCfg.nSubsteps = 200;
% T=100*0.2=20; dt=0.2/200=0.001.

convergenceCriterion = "norm(Wnew-Wold,'fro') <= relTol*norm(Wold,'fro')";
simulationMethod = "Euler-Maruyama; one scalar Wiener process per path";
matlabVersion = version;
dependencyFiles = {fullfile('dependencies', 'riccati.m'); ...
    fullfile('dependencies', 'sample_path.m')};

%% Solve, simulate, and save both cases
pathProbes = cell(2, 1);
pathResults = cell(2, 1);
pathFigures = gobjects(2, 1);
pathAxes = gobjects(2, 1);
runTimer = tic;
for bIndex = 1:2
    pathModel = pathModels{bIndex};
    fprintf('\n=== %s ===\n', pathCaseNames{bIndex});
    pathProbe = riccati(pathModel, solverCfg);
    pathProbes{bIndex} = pathProbe;
    if ~pathProbe.converged
        error('sample_paths:RiccatiNotConverged', ...
            '%s: %s.', pathCaseNames{bIndex}, pathProbe.status);
    end
    [pathFigures(bIndex), pathAxes(bIndex), pathResults{bIndex}] = ...
        sample_path(pathModel, pathProbe.Wfinal, pathCfg);
    set(pathFigures(bIndex), 'Name', pathCaseNames{bIndex}, ...
        'NumberTitle', 'off');
    savefig(pathFigures(bIndex), ...
        fullfile(outputDir, [pathCaseNames{bIndex} '.fig']));
    exportgraphics(pathAxes(bIndex), ...
        fullfile(outputDir, [pathCaseNames{bIndex} '.png']), 'Resolution', 300);
end
runSeconds = toc(runTimer);
save(fullfile(outputDir, 'results.mat'), ...
    'pathModels', 'solverCfg', 'pathCfg', 'pathCaseNames', ...
    'pathProbes', 'pathResults', 'convergenceCriterion', ...
    'simulationMethod', 'matlabVersion', 'dependencyFiles', 'runSeconds');
fprintf('\nSample paths saved to:\n%s\n', outputDir);
