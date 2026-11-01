function comparison = run_A3_relTol_comparison
%RUN_A3_RELTOL_COMPARISON Rerun A_3 boundary curves for four stopping tolerances.
% Uses the saved comparison settings; each candidate starts from P=0.

    outputDir = fileparts(mfilename('fullpath'));
    projectDir = fileparts(fileparts(outputDir));
    addpath(fullfile(projectDir, 'dependencies'), '-begin');
    reference = load(fullfile(outputDir, 'results.mat'));

    aIndex = reference.aIndex;
    A3 = reference.A3;
    Bset = reference.Bset;
    Rset = reference.Rset;
    L = reference.L;
    P = zeros(size(A3,1));
    tauSet = reference.tauSet;
    baseSolverCfg = reference.solverConfigs{1};
    if isfield(reference,'searchConfigs')
        searchCfg = reference.searchConfigs{1};
    else
        searchCfg = reference.searchCfg;
    end
    convergenceCriterion = reference.convergenceCriterion;
    relTols = [1e-3, 1e-4, 1e-5, 1e-6];
    matlabVersion = version;

    nTols = numel(relTols);
    solverConfigs = cell(nTols,1);
    searchConfigs = repmat({searchCfg},nTols,1);
    sigmaCurveI = nan(numel(tauSet{1}),nTols);
    sigmaCurveOnes = nan(numel(tauSet{2}),nTols);
    curveBoundaryI = cell(size(sigmaCurveI));
    curveBoundaryOnes = cell(size(sigmaCurveOnes));
    scanSeconds = zeros(nTols,2);
    model = struct('A',A3,'B',Bset{1},'L',L,'R',Rset{1}, ...
        'P',P,'tau',0,'sigma',0);
    timer = tic;
    fprintf('A_3 comparison: four tolerances, %d boundary jobs.\n', ...
        nTols*(numel(tauSet{1})+numel(tauSet{2})));
    fprintf('The scan labels its single supplied matrix A_1; this matrix is A_3.\n');

    for k=1:nTols
        cfg = baseSolverCfg;
        cfg.relTol = relTols(k);
        solverConfigs{k} = cfg;
        fprintf('\n=== A_3 | B = I_5 | epsilon = %.1e ===\n',cfg.relTol);
        model.B = Bset{1};
        model.R = Rset{1};
        [sigmaCurveI(:,k),curveBoundaryI(:,k),scanSeconds(k,1)] = ...
            scan_sigma_curves({A3},tauSet{1},model,cfg,searchCfg);
        fprintf('\n=== A_3 | B = ones(5,1) | epsilon = %.1e ===\n',cfg.relTol);
        model.B = Bset{2};
        model.R = Rset{2};
        [sigmaCurveOnes(:,k),curveBoundaryOnes(:,k),scanSeconds(k,2)] = ...
            scan_sigma_curves({A3},tauSet{2},model,cfg,searchCfg);
    end

    sigmaSquaredI = sigmaCurveI.^2;
    sigmaSquaredOnes = sigmaCurveOnes.^2;
    relativeDifferenceI = (sigmaSquaredI-sigmaSquaredI(:,end))./sigmaSquaredI(:,end);
    relativeDifferenceOnes = ...
        (sigmaSquaredOnes-sigmaSquaredOnes(:,end))./sigmaSquaredOnes(:,end);
    comparisonSummary = struct('referenceRelTol',relTols(end), ...
        'nDifferentI',sum(sigmaSquaredI~=sigmaSquaredI(:,end),1), ...
        'nDifferentOnes',sum(sigmaSquaredOnes~=sigmaSquaredOnes(:,end),1), ...
        'maxAbsRelativeDifferenceI',max(abs(relativeDifferenceI),[],1), ...
        'maxAbsRelativeDifferenceOnes',max(abs(relativeDifferenceOnes),[],1));
    runSeconds = toc(timer);
    save(fullfile(outputDir,'results.mat'), ...
        'aIndex','A3','Bset','Rset','L','P','tauSet','relTols', ...
        'solverConfigs','searchConfigs','convergenceCriterion', ...
        'sigmaCurveI','sigmaCurveOnes','sigmaSquaredI','sigmaSquaredOnes', ...
        'curveBoundaryI','curveBoundaryOnes','scanSeconds','runSeconds', ...
        'relativeDifferenceI','relativeDifferenceOnes','comparisonSummary','matlabVersion');

    % The adjacent plotting script keeps figure styling separate from the scans.
    run(fullfile(outputDir,'plot_A3_relTol_comparison.m'));
    comparison = load(fullfile(outputDir,'results.mat'));
end
