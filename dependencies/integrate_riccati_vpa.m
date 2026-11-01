function result = integrate_riccati_vpa(model, cfg, rhsFunction)
%INTEGRATE_RICCATI_VPA Integrate the continuous Riccati equation with VPA.
% Uses adaptive Dormand-Prince 5(4) steps. The caller sets digits and passes
% a high-precision model; rhsFunction returns [dW, failure, residualScale].
% odeRelTol/odeAbsTol control integration error; relTol controls stationarity.
% odeMaxStep caps the step size and maxIter caps the accepted step count.

    timer = tic;
    W = (model.P + model.P.') / 2;
    time = vpa(0);
    nUpdates = 0;
    status = 'inconclusive_max_iter';
    toVpa = @(value) vpa(sym(value, 'f'), cfg.precisionDigits);
    maxTime = toVpa(cfg.odeMaxTime);
    maxStep = toVpa(cfg.odeMaxStep);
    absTol = toVpa(cfg.odeAbsTol);
    relTol = toVpa(cfg.odeRelTol);
    normLimit = toVpa(cfg.maxWNorm);
    stoppingTolerance = toVpa(cfg.relTol);
    step = min(maxStep, maxTime);
    [stageWeights, errorWeights] = dormand_prince_weights;

    hasQueue = isfield(cfg, 'progressQueue');
    if isfield(cfg, 'progressInterval')
        reportInterval = cfg.progressInterval;
    else
        reportInterval = 1;
    end
    lastReportTime = -Inf;

    [dW, failure, residualScale] = rhsFunction(W, model);
    if ~isempty(failure)
        normalizedResidual = NaN;
        status = failure;
        result = finish_result;
        return;
    end
    normalizedResidual = residual_ratio(norm(dW, 'fro'), residualScale);
    if normalizedResidual <= stoppingTolerance
        status = 'converged';
        result = finish_result;
        return;
    end
    report_progress;

    while time < maxTime && nUpdates < cfg.maxIter
        remainingTime = maxTime - time;
        step = min(step, min(maxStep, remainingTime));
        isFinalStep = isAlways(step == remainingTime, 'Unknown', 'error');
        smallestStep = vpa(10)^(6 - cfg.precisionDigits) * max(1, abs(time));
        if step <= 0 || (step < smallestStep && ~isFinalStep)
            status = 'inconclusive_ode_step_underflow';
            break;
        end

        [candidate, candidateDW, candidateScale, errorEstimate, failure] = ...
            rk45_step(W, dW, step, model, rhsFunction, ...
                stageWeights, errorWeights);
        if ~isempty(failure)
            % A rejected intermediate stage is not an accepted trajectory value.
            step = step / 2;
            report_progress;
            continue;
        end

        errorScale = absTol + relTol * max(abs(W), abs(candidate));
        errorRatio = max(abs(errorEstimate(:)) ./ errorScale(:));
        if ~isfinite(errorRatio)
            step = step / 2;
            report_progress;
            continue;
        end

        if errorRatio <= 1
            if norm(candidate, 'fro') >= normLimit
                % Locate the norm limit by recomputing shorter RK steps.
                [W, fraction] = norm_limit_point(W, dW, candidate, step, ...
                    model, rhsFunction, stageWeights, errorWeights, normLimit);
                time = time + fraction * step;
                nUpdates = nUpdates + 1;
                [dW, failure, residualScale] = rhsFunction(W, model);
                if isempty(failure)
                    normalizedResidual = residual_ratio(norm(dW, 'fro'), residualScale);
                else
                    normalizedResidual = NaN;
                end
                status = 'inconclusive_W_norm_limit';
                break;
            end

            W = candidate;
            dW = candidateDW;
            time = time + step;
            if isFinalStep
                time = maxTime;
            end
            nUpdates = nUpdates + 1;
            normalizedResidual = residual_ratio(norm(dW, 'fro'), candidateScale);
            if normalizedResidual <= stoppingTolerance
                status = 'converged';
                break;
            end
        end

        if errorRatio == 0
            factor = vpa(5);
        else
            factor = (sym(9) / 10) * errorRatio^(-sym(1) / 5);
            factor = min(5, max(sym(1) / 5, factor));
        end
        if errorRatio > 1
            factor = min(sym(1) / 2, factor);
        end
        step = min(maxStep, step * factor);
        report_progress;
    end

    if strcmp(status, 'inconclusive_max_iter') && time >= maxTime
        status = 'inconclusive_max_time';
    end
    result = finish_result;

    function report_progress
        elapsed = toc(timer);
        if elapsed - lastReportTime < reportInterval
            return;
        end
        lastReportTime = elapsed;
        stopRatio = double(normalizedResidual / stoppingTolerance);
        if hasQueue
            send(cfg.progressQueue, struct('event', 'ode', ...
                'jobIndex', cfg.progressJobIndex, ...
                'sigmaSquared', double(model.sigma^2), ...
                'step', nUpdates, 'stopRatio', stopRatio, ...
                'odeTime', double(time)));
        else
            fprintf(['ODE t %.6g/%.6g | steps %d | normalizedResidual %.3e' ...
                ' (tol %.1e) | elapsed %.1fs\n'], ...
                double(time), cfg.odeMaxTime, nUpdates, ...
                double(normalizedResidual), cfg.relTol, elapsed);
        end
    end

    function packed = finish_result
        packed.Wfinal = W;
        packed.converged = strcmp(status, 'converged');
        packed.status = status;
        packed.nUpdates = nUpdates;
        packed.normalizedResidual = normalizedResidual;
        packed.method = 'continuous_riccati';
        packed.finalTime = double(time);
        packed.iterationSeconds = toc(timer);
        if ~hasQueue
            fprintf(['%s | steps %d | checkedResidual %.3e | finalTime %.6g' ...
                ' | integration time %.3fs\n'], ...
                status, nUpdates, double(normalizedResidual), ...
                packed.finalTime, packed.iterationSeconds);
        end
    end
end

function [candidate, candidateDW, candidateScale, errorEstimate, failure] = ...
        rk45_step(W, dW, step, model, rhsFunction, a, errorWeights)
    n = size(W, 1);
    derivatives = sym(zeros(numel(W), 7));
    derivatives(:, 1) = dW(:);
    candidate = [];
    candidateDW = [];
    candidateScale = NaN;
    errorEstimate = [];
    failure = '';
    for stage = 2:7
        stageVector = W(:) + step * ...
            (derivatives(:, 1:stage-1) * a(stage, 1:stage-1).');
        stageW = reshape(stageVector, n, n);
        stageW = (stageW + stageW.') / 2;
        [stageDW, failure, stageScale] = rhsFunction(stageW, model);
        if ~isempty(failure)
            return;
        end
        derivatives(:, stage) = stageDW(:);
    end
    % The seventh stage is the fifth-order candidate and its fresh derivative.
    candidate = stageW;
    candidateDW = stageDW;
    candidateScale = stageScale;
    errorEstimate = reshape(step * (derivatives * errorWeights), n, n);
end

function [a, errorWeights] = dormand_prince_weights
    % Exact rational coefficients, converted at the caller's working precision.
    a = sym(zeros(7));
    a(2, 1) = sym(1) / 5;
    a(3, 1:2) = sym([3, 9]) / 40;
    a(4, 1:3) = sym([44, -56, 32]) ./ sym([45, 15, 9]);
    a(5, 1:4) = sym([19372, -25360, 64448, -212]) ./ ...
        sym([6561, 2187, 6561, 729]);
    a(6, 1:5) = sym([9017, -355, 46732, 49, -5103]) ./ ...
        sym([3168, 33, 5247, 176, 18656]);
    a(7, 1:6) = sym([35, 0, 500, 125, -2187, 11]) ./ ...
        sym([384, 1, 1113, 192, 6784, 84]);
    fifthOrder = a(7, :).';
    fourthOrder = sym([5179; 0; 7571; 393; -92097; 187; 1]) ./ ...
        sym([57600; 1; 16695; 640; 339200; 2100; 40]);
    errorWeights = vpa(fifthOrder - fourthOrder);
    a = vpa(a);
end

function ratio = residual_ratio(residual, scale)
    if scale > 0
        ratio = residual / scale;
    elseif residual == 0
        ratio = vpa(0);
    else
        ratio = Inf;
    end
end

function [crossing, fraction] = norm_limit_point(previous, previousDW, ...
        candidate, step, model, rhsFunction, a, errorWeights, limit)
    lower = vpa(0);
    upper = vpa(1);
    crossing = candidate;
    for iteration = 1:50
        midpoint = (lower + upper) / 2;
        [midpointW, ~, ~, ~, failure] = rk45_step(previous, previousDW, ...
            midpoint * step, model, rhsFunction, a, errorWeights);
        if ~isempty(failure)
            % Keep the valid upper state if a refined trial cannot be evaluated.
            break;
        end
        if norm(midpointW, 'fro') >= limit
            upper = midpoint;
            crossing = midpointW;
        else
            lower = midpoint;
        end
    end
    fraction = upper;
    crossing = (crossing + crossing.') / 2;
end
