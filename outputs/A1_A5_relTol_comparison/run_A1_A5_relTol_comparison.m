function comparison = run_A1_A5_relTol_comparison
%RUN_A1_A5_RELTOL_COMPARISON Rerun four tolerances for A_1 through A_5.
% Load the adjacent saved settings and bisect sigma^2 in [0,4] for every case.
% Requires MATLAB and Parallel Computing Toolbox. Each solve starts at P=0.

    outputDir = fileparts(mfilename('fullpath'));
    projectDir = fileparts(fileparts(outputDir));
    addpath(fullfile(projectDir, 'dependencies'), '-begin');
    reference = load(fullfile(outputDir, 'results.mat'));

    Aset = reference.Aset;
    Bset = reference.Bset;
    Rset = reference.Rset;
    L = reference.L;
    P = reference.P;
    tauSet = reference.tauSet;
    relTols = [1e-3, 1e-4, 1e-5, 1e-6];
    searchCfg = reference.searchCfg;
    solverConfigs = repmat(reference.solverConfigs(1), size(relTols));
    convergenceCriterion = reference.convergenceCriterion;
    searchVariable = 'sigmaSquared';
    matlabVersion = version;
    timer = tic;

    for k = 1:numel(relTols)
        solverConfigs{k}.relTol = relTols(k);
    end
    fprintf('Five-matrix comparison: %d boundary jobs, bisection in sigma^2.\n', ...
        numel(Aset) * sum(cellfun(@numel, tauSet)) * numel(relTols));

    for bIndex = 1:numel(Bset)
        tau = tauSet{bIndex};
        model = struct('A', Aset{1}, 'B', Bset{bIndex}, 'L', L, ...
            'R', Rset{bIndex}, 'P', P, 'tau', 0, 'sigma', 0);
        sigmaCurve = nan(numel(tau), numel(Aset), numel(relTols));
        boundaries = cell(size(sigmaCurve));
        scanSeconds = zeros(size(relTols));

        for k = 1:numel(relTols)
            fprintf('\n=== Input case %d | epsilon = %.1e ===\n', bIndex, relTols(k));
            [sigmaCurve(:,:,k), boundaries(:,:,k), scanSeconds(k)] = ...
                scan_sigma_curves(Aset, tau, model, solverConfigs{k}, searchCfg);
        end

        sigmaSquared = cellfun(@(boundary) boundary.sigmaSquareLow, boundaries);
        referenceValues = sigmaSquared(:,:,end);
        comparable = all(isfinite(sigmaSquared) & sigmaSquared > 0, 3);
        relativeDifference = (sigmaSquared - referenceValues) ./ referenceValues;
        maxRelativeDifference = max(abs(relativeDifference), [], 3);
        maxRelativeDifference(~comparable) = NaN;
        different = comparable & any(sigmaSquared ~= referenceValues, 3);
        [tauIndex, aIndex] = find(different);
        values = reshape(sigmaSquared, [], numel(relTols));
        changedRows = table(aIndex, tau(tauIndex), values(different(:),:), ...
            100 * maxRelativeDifference(different), ...
            'VariableNames', {'AIndex', 'Tau', 'SigmaSquared', 'MaxRelativeDifferencePercent'});

        cases(bIndex) = struct('sigmaCurve', sigmaCurve, 'sigmaSquared', sigmaSquared, ...
            'boundaries', {boundaries}, 'scanSeconds', scanSeconds, ...
            'comparable', comparable, 'changedRows', changedRows);
    end

    runSeconds = toc(timer);
    resultFile = fullfile(outputDir, 'results.mat');
    save(resultFile, 'Aset', 'Bset', 'Rset', 'L', 'P', 'tauSet', 'relTols', ...
        'solverConfigs', 'searchCfg', 'searchVariable', 'convergenceCriterion', ...
        'cases', 'matlabVersion', 'runSeconds');
    run(fullfile(outputDir, 'plot_A1_A5_relTol_comparison.m'));
    comparison = load(resultFile);
end
