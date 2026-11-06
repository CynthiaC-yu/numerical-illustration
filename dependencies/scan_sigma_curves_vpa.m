function [sigmaCurve, curveBoundary, curveSeconds] = ...
    scan_sigma_curves_vpa(Aset, tauValues, model, solverCfg, searchCfg, checkpointDir)
%SCAN_SIGMA_CURVES_VPA Run independent VPA boundary searches with parfor.
% Rows correspond to tauValues; columns correspond to Aset. Values are sigma.
% An optional checkpoint directory stores each completed boundary as a MAT file.

    if nargin < 6
        checkpointDir = '';
    end
    if ~isempty(checkpointDir) && ~isfolder(checkpointDir)
        mkdir(checkpointDir);
    end
    tauTest = tauValues(:);
    nTau = numel(tauTest);
    nA = numel(Aset);
    modelJobs = cell(nTau, nA);
    for aIndex = 1:nA
        for tauIndex = 1:nTau
            modelJob = model;
            modelJob.A = Aset{aIndex};
            modelJob.tau = tauTest(tauIndex);
            modelJobs{tauIndex, aIndex} = modelJob;
        end
    end
    nJobs = numel(modelJobs);
    curveBoundary = cell(size(modelJobs));
    sigmaCurve = nan(size(modelJobs));

    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('Processes', 'IdleTimeout', Inf);
    end
    if ~isa(pool, 'parallel.ProcessPool')
        error('scan_sigma_curves_vpa:PoolType', 'VPA scans require a Processes pool.');
    end
    fprintf('\nParallel VPA scan: %d jobs on %d workers.\n', nJobs, pool.NumWorkers);

    if ~isfield(solverCfg, 'progressInterval')
        solverCfg.progressInterval = 5;
    end
    % Keep saved configuration free of the queue and job-specific metadata.
    savedSolverCfg = solverCfg;
    progressQueue = parallel.pool.DataQueue;
    solverCfg.progressQueue = progressQueue;
    startedJobs = false(nJobs, 1);
    finishedJobs = false(nJobs, 1);
    probeCounts = zeros(nJobs, 1);
    lastPrintTime = -Inf;
    curveTimer = tic;
    progressListener = afterEach(progressQueue, @display_progress);
    listenerCleanup = onCleanup(@() delete(progressListener));

    parfor jobIndex = 1:nJobs
        modelJob = modelJobs{jobIndex};
        jobCfg = solverCfg;
        jobCfg.progressJobIndex = jobIndex;
        send(progressQueue, struct('event', 'started', 'jobIndex', jobIndex));
        pointBoundary = find_sigma_boundary_vpa(modelJob, jobCfg, searchCfg);

        if ~isempty(checkpointDir)
            [tauIndex, aIndex] = ind2sub([nTau, nA], jobIndex);
            save_boundary(fullfile(checkpointDir, sprintf('job_%03d.mat', jobIndex)), ...
                pointBoundary, modelJob, savedSolverCfg, searchCfg, aIndex, tauIndex);
        end
        curveBoundary{jobIndex} = pointBoundary;
        sigmaCurve(jobIndex) = sqrt(pointBoundary.sigmaSquareLow);
        send(progressQueue, struct('event', 'finished', 'jobIndex', jobIndex, ...
            'sigmaSquared', pointBoundary.sigmaSquareLow));
    end

    drawnow;
    while progressQueue.QueueLength > 0
        drawnow;
    end
    curveSeconds = toc(curveTimer);
    fprintf('Parallel VPA scan finished: %d/%d jobs in %.2fs.\n', ...
        nJobs, nJobs, curveSeconds);
    nInvalid = nnz(~isfinite(sigmaCurve));
    if nInvalid > 0
        fprintf('%d jobs have no valid search bracket; inspect curveBoundary.\n', nInvalid);
    end

    function display_progress(message)
        index = message.jobIndex;
        switch message.event
            case 'started'
                startedJobs(index) = true;
            case 'probe'
                probeCounts(index) = probeCounts(index) + 1;
            case 'finished'
                finishedJobs(index) = true;
        end
        elapsed = toc(curveTimer);
        nDone = nnz(finishedJobs);
        if elapsed - lastPrintTime < solverCfg.progressInterval && nDone < nJobs
            return;
        end
        lastPrintTime = elapsed;
        nActive = nnz(startedJobs & ~finishedJobs);
        if nDone > 0
            etaText = ['~' format_duration(elapsed * (nJobs - nDone) / nDone)];
        else
            etaText = '--';
        end
        fprintf('[%d/%d | %.1f%%] active %d | elapsed %s | ETA %s\n', ...
            nDone, nJobs, 100 * nDone / nJobs, nActive, ...
            format_duration(elapsed), etaText);

        [tauIndexNow, aIndexNow] = ind2sub([nTau, nA], index);
        fprintf('  A_%d | tau = %.3g', aIndexNow, tauTest(tauIndexNow));
        switch message.event
            case 'started'
                fprintf(' | starting\n');
            case 'probe'
                fprintf(' | probe %d | sigma^2 = %.6g\n', ...
                    probeCounts(index), message.sigmaSquared);
            case 'iteration'
                fprintf([' | probe %d | sigma^2 = %.6g | iter %d/%d' ...
                    ' | stopRatio %.3e (pass <= 1)\n'], ...
                    probeCounts(index), message.sigmaSquared, ...
                    message.step, solverCfg.maxIter, message.stopRatio);
            case 'ode'
                fprintf([' | probe %d | sigma^2 = %.6g | ODE t = %.4g/%.4g' ...
                    ' | stopRatio %.3e (pass <= 1)\n'], ...
                    probeCounts(index), message.sigmaSquared, ...
                    message.odeTime, solverCfg.odeMaxTime, message.stopRatio);
            case 'finished'
                if isfinite(message.sigmaSquared)
                    fprintf(' | done | sigma^2 = %.6g\n', message.sigmaSquared);
                else
                    fprintf(' | done | no valid search bracket\n');
                end
        end
    end
end

function save_boundary(fileName, boundary, model, solverCfg, searchCfg, aIndex, tauIndex)
    save(fileName, 'boundary', 'model', 'solverCfg', 'searchCfg', 'aIndex', 'tauIndex');
end

function text = format_duration(seconds)
    if seconds < 60
        text = sprintf('%.1fs', seconds);
    elseif seconds < 3600
        text = sprintf('%.1fmin', seconds / 60);
    else
        text = sprintf('%.2fh', seconds / 3600);
    end
end
