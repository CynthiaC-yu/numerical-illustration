%% Numerical illustration with the min criterion and VPA arithmetic
clear;
clc;

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'dependencies'), '-begin');

%% 1. Define the models and sampling grids
n = 5;
Aset = cell(5, 1);
for i = 1:5
    Aset{i} = ((i-1)/5) * eye(n) + diag((0:n-1)/10);
end
Bset = {eye(n); ones(n, 1)};
Rset = {eye(n); 1};
L = eye(n);
tauSet = {(0:0.2:10).'; (0:0.2:6).'};

modelOnes = struct('A', Aset{1}, 'B', Bset{2}, ...
    'L', L, 'R', Rset{2}, 'P', zeros(n), 'tau', 0, 'sigma', 0);
modelI = modelOnes;
modelI.B = Bset{1};
modelI.R = Rset{1};

%% 2. Configure the solver and boundary search
solverCfg.precisionDigits = 50;
solverCfg.Ns = 201;
solverCfg.maxIter = 2e5;
solverCfg.maxWNorm = 1e12;
% norm(Wnew-Wold,'fro') <= relTol * min(1,norm(Wold,'fro')).
solverCfg.relTol = 1e-3;
solverCfg.progressInterval = 5;

% At tau=0, integrate in VPA and check the min-scaled derivative residual.
solverCfg.odeSolver = 'rk45_vpa';
solverCfg.odeMaxTime = 5000;
solverCfg.odeMaxStep = 0.1;
solverCfg.odeRelTol = 1e-9;
solverCfg.odeAbsTol = 1e-11;
searchCfg.sigmaSquaredRelTol = 0.02;

%% 3. Start or reuse a process pool
% Set nWorkers=24 on the workstation; [] uses MATLAB's default pool size.
nWorkers = [];
pool = gcp('nocreate');
if isempty(pool)
    cluster = parcluster('Processes');
    if isempty(nWorkers)
        pool = parpool(cluster, 'IdleTimeout', Inf);
    else
        cluster.NumWorkers = nWorkers;
        pool = parpool(cluster, nWorkers, 'IdleTimeout', Inf);
    end
elseif ~isa(pool, 'parallel.ProcessPool')
    error('main_vpa:PoolType', 'VPA scans require a Processes pool.');
elseif ~isempty(nWorkers) && pool.NumWorkers ~= nWorkers
    error('main_vpa:PoolSize', ...
        'The existing pool has %d workers; requested %d.', pool.NumWorkers, nWorkers);
end
actualWorkers = pool.NumWorkers;

%% 4. Save settings before the long scan
runStamp = char(datetime('now', 'TimeZone', 'America/Chicago', ...
    'Format', 'yyyyMMdd_HHmmss_SSS'));
outputDir = fullfile(projectDir, 'outputs', [runStamp '_vpa']);
mkdir(outputDir);
save(fullfile(outputDir, 'settings.mat'), ...
    'Aset', 'Bset', 'Rset', 'L', 'tauSet', 'modelI', 'modelOnes', ...
    'solverCfg', 'searchCfg', 'nWorkers', 'actualWorkers');
fprintf('\nVPA run: %d digits, %d workers, %d boundary jobs.\n', ...
    solverCfg.precisionDigits, actualWorkers, ...
    numel(Aset) * (numel(tauSet{1}) + numel(tauSet{2})));
fprintf('Output folder:\n%s\n', outputDir);

%% 5. Scan B = ones(5,1) first; save each completed boundary
fprintf('\n=== B = ones(5,1) ===\n');
[sigmaCurve, curveBoundary, curveSeconds] = scan_sigma_curves_vpa( ...
    Aset, tauSet{2}, modelOnes, solverCfg, searchCfg, ...
    fullfile(outputDir, 'checkpoints', 'B_ones5'));
save(fullfile(outputDir, 'B_ones5_results.mat'), ...
    'sigmaCurve', 'curveBoundary', 'curveSeconds');

fprintf('\n=== B = I_5 ===\n');
[sigmaCurveI, curveBoundaryI, curveSecondsI] = scan_sigma_curves_vpa( ...
    Aset, tauSet{1}, modelI, solverCfg, searchCfg, ...
    fullfile(outputDir, 'checkpoints', 'B_I5'));
save(fullfile(outputDir, 'B_I5_results.mat'), ...
    'sigmaCurveI', 'curveBoundaryI', 'curveSecondsI');

save(fullfile(outputDir, 'results.mat'), ...
    'Aset', 'Bset', 'Rset', 'L', 'tauSet', 'solverCfg', 'searchCfg', ...
    'actualWorkers', 'sigmaCurveI', 'curveBoundaryI', 'curveSecondsI', ...
    'sigmaCurve', 'curveBoundary', 'curveSeconds');

%% 6. Plot and export both figures
[figI, axI] = plot_sigma_curves(tauSet{1}, sigmaCurveI, '', 2);
ylim(axI, [1e-4, 2]);
[figOnes, axOnes] = plot_sigma_curves(tauSet{2}, sigmaCurve, '', 1);
savefig(figI, fullfile(outputDir, 'B_I5.fig'));
exportgraphics(figI, fullfile(outputDir, 'B_I5.png'), 'Resolution', 300);
savefig(figOnes, fullfile(outputDir, 'B_ones5.fig'));
exportgraphics(figOnes, fullfile(outputDir, 'B_ones5.png'), 'Resolution', 300);
fprintf('\nResults and figures saved to:\n%s\n', outputDir);
