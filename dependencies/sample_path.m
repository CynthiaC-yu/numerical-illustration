function [fig, ax, pathResult] = sample_path(model, W, pathCfg)
%SAMPLE_PATH Simulate and plot squared state norms for one chosen model.
% The caller supplies W from a converged Riccati solve for this same model.
% pathCfg requires x0, rngSeed, nPaths, nShownPaths, nSampleIntervals,
% and nSubsteps. Each call creates one figure and returns the trajectory data.
% Model selection, Riccati iteration, and file saving belong to the caller.

    A = model.A;
    B = model.B;
    L = model.L;
    R = model.R;
    tau = model.tau;
    sigma = model.sigma;
    n = size(A, 1);
    m = size(B, 2);

    validateattributes(tau, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
    validateattributes(W, {'numeric'}, {'real', 'finite', 'size', [n n]});
    validateattributes(pathCfg.x0, {'numeric'}, {'real', 'finite', 'size', [n 1]});
    countFields = {'nPaths', 'nShownPaths', 'nSampleIntervals', 'nSubsteps'};
    for k = 1:numel(countFields)
        validateattributes(pathCfg.(countFields{k}), {'numeric'}, ...
            {'real', 'finite', 'scalar', 'integer', 'positive'});
    end
    if pathCfg.nShownPaths > pathCfg.nPaths
        error('sample_path:PathCount', 'nShownPaths must not exceed nPaths.');
    end

    %% Compute the gain K(s) from the supplied W
    nSubsteps = pathCfg.nSubsteps;
    sGrid = linspace(0, tau, nSubsteps + 1);
    Kgrid = zeros(m, n, numel(sGrid));
    for j = 1:numel(sGrid)
        s = sGrid(j);
        expARemainingB = expm(A * (tau - s)) * B;
        Q = R + sigma^2 * expARemainingB.' * W * expARemainingB;
        Q = (Q + Q.') / 2;
        [Qfactor, cholFlag] = chol(Q, 'lower');
        if cholFlag ~= 0
            error('sample_path:NonpositiveQ', ...
                'Q(tau-s) must be positive definite when computing K(s).');
        end
        rhs = B.' * expm(-A.' * s) * (W - tau * L);
        Kgrid(:, :, j) = -(Qfactor.' \ (Qfactor \ rhs));
    end

    % Euler-Maruyama uses the left endpoint of each substep.
    BKgrid = zeros(n, n, nSubsteps);
    for j = 1:nSubsteps
        BKgrid(:, :, j) = B * Kgrid(:, :, j);
    end

    %% Simulate the sample paths
    previousRng = rng;
    rngCleanup = onCleanup(@() rng(previousRng));
    rng(pathCfg.rngSeed, 'twister');

    nPaths = pathCfg.nPaths;
    nShownPaths = pathCfg.nShownPaths;
    nSteps = pathCfg.nSampleIntervals * nSubsteps;
    dt = tau / nSubsteps;
    tGrid = (0:nSteps).' * dt;
    X = repmat(pathCfg.x0, 1, nPaths);

    meanNormSquared = zeros(nSteps + 1, 1);
    shownNormSquared = zeros(nSteps + 1, nShownPaths);
    firstStatePath = zeros(nSteps + 1, n);
    initialNormSquared = sum(X.^2, 1);
    meanNormSquared(1) = mean(initialNormSquared);
    shownNormSquared(1, :) = initialNormSquared(1:nShownPaths);
    firstStatePath(1, :) = X(:, 1).';

    globalStep = 0;
    for intervalIndex = 1:pathCfg.nSampleIntervals
        % Hold x_{k*tau} fixed throughout this sampling interval.
        Xsample = X;
        for j = 1:nSubsteps
            globalStep = globalStep + 1;
            controlEffect = BKgrid(:, :, j) * Xsample;
            % One scalar Brownian increment for each independent path.
            dW = sqrt(dt) * randn(1, nPaths);
            X = X + (A * X + controlEffect) * dt ...
                + sigma * (controlEffect .* dW);

            currentNormSquared = sum(X.^2, 1);
            meanNormSquared(globalStep + 1) = mean(currentNormSquared);
            shownNormSquared(globalStep + 1, :) = ...
                currentNormSquared(1:nShownPaths);
            firstStatePath(globalStep + 1, :) = X(:, 1).';
        end
        if any(~isfinite(X(:)))
            error('sample_path:NonfiniteState', ...
                'The simulated state became nonfinite. Check the model, gain, and time step.');
        end
    end
    samplingIndices = 1:nSubsteps:(nSteps + 1);

    %% Return the model, settings, gains, and trajectory data
    pathResult.model = model;
    pathResult.W = W;
    pathResult.config = pathCfg;
    pathResult.dt = dt;
    pathResult.sGrid = sGrid;
    pathResult.Kgrid = Kgrid;
    pathResult.tGrid = tGrid;
    pathResult.meanNormSquared = meanNormSquared;
    pathResult.shownNormSquared = shownNormSquared;
    pathResult.firstStatePath = firstStatePath;
    pathResult.finalStates = X;
    pathResult.samplingIndices = samplingIndices;

    %% Plot the individual squared norms and their Monte Carlo mean
    fig = figure('Color', 'w', 'Units', 'pixels', ...
        'Position', [100 100 968 660]);
    ax = axes('Parent', fig);
    sampleHandles = semilogy(ax, tGrid, max(shownNormSquared, realmin), ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    hold(ax, 'on');
    meanHandle = semilogy(ax, tGrid, max(meanNormSquared, realmin), ...
        'Color', [0 0.447 0.741], 'LineWidth', 2.2);
    samplingHandle = semilogy(ax, tGrid(samplingIndices), ...
        max(meanNormSquared(samplingIndices), realmin), 'o', ...
        'Color', [0.85 0.325 0.098], 'MarkerSize', 4, ...
        'MarkerFaceColor', [0.85 0.325 0.098]);
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Helvetica', 'FontSize', 20, ...
        'LineWidth', 1.1, 'TickDir', 'in');
    xlabel(ax, 't', 'FontSize', 22);
    ylabel(ax, '||x_t||_2^2', 'Interpreter', 'tex', 'FontSize', 22);
    legend(ax, [sampleHandles(1), meanHandle, samplingHandle], ...
        {'Individual paths', 'Monte Carlo mean', 'Sampling instants'}, ...
        'Location', 'best', 'FontSize', 18);
end
