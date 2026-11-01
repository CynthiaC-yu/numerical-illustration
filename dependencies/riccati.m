function result = riccati(model, cfg)
    % ----- Check if riccati solver received a complete signature.
    if nargin < 2
        error('riccati:MissingInput', ...
            'Both model and cfg must be provided.')
    end

    % ----- Check if model and config input is complete and valid.
    validate_model(model, cfg.maxWNorm);
    validate_config(cfg, model.tau);

    % ----- Loading model, config, and initialize variables
    A = model.A;
    B = model.B;
    L = model.L;
    R = model.R;
    tau = model.tau;
    sigmaSquared = model.sigma^2;
    P = model.P;
    % Initialize dimension parameters
    n = size(A, 1);
    m = size(B, 2);
    % Initialize W (We have checked P is symmtric within numerical tolerance)
    W = 0.5 * P + 0.5 * P.';

    % During a parallel scan, report through the client's progress queue.
    hasProgressQueue = isfield(cfg, 'progressQueue');
    if hasProgressQueue
        report_progress(cfg, 'probe', model, 0, NaN, NaN);
    end

    % ----- For Tau > 0 -----
    if tau > 0
        % ----- Precomputing section ------
        % Creates grid and quadrature weights for numerical integral.
        Ns = cfg.Ns;
        sGrid = linspace(0, tau, Ns);
        h = tau / (Ns - 1);
        quadratureWeights = h * ones(1, Ns);
        quadratureWeights(1) = h / 2;
        quadratureWeights(end) = h / 2;
        % Compute expAsB Grid (i.e. exp(As)B)
        expAsBGrid = zeros(n, m, Ns);
        for j = 1:Ns
            expAsBGrid(:, :, j) = expm(A * sGrid(j)) * B;
        end
        % Compute expATau which used by every Riccati update (i.e. f(W)).
        expATau = expm(A * tau);
        % Assigning vars
        precomputed.expAsBGrid = expAsBGrid;
        precomputed.quadratureWeights = quadratureWeights;
        precomputed.expATau = expATau;

        % ----- The loop for updating W -----
        % Initiate key variables
        status = 'inconclusive_max_iter';
        nUpdates = 0;
        normalizedResidual = NaN;
        updateResidual = NaN;
        relativeUpdate = NaN;
        % Inititate progress recorder
        progressTimer = tic;
        lastProgressTime = 0;
        if hasProgressQueue
            progressInterval = cfg.progressInterval;
        else
            progressInterval = 1;
        end
        for updateIndex = 1:cfg.maxIter
            % WNew = f(W) and then check if there's any failure before
            % proceeding
            [Wnew, failure] = compute_f(W, model, precomputed, cfg.maxWNorm);
            if ~isempty(failure)
                status = failure;
                break;
            end
            %  Compute |f^{N}(P) - f^{N-1}(P)|_F
            updateResidual = norm(Wnew - W, 'fro');
            oldWNorm = norm(W, 'fro');
            % Normalize the difference: |f^{N}(P) - f^{N-1}(P)|_F / |f^{N-1}(P)|_F
            if oldWNorm > 0
                relativeUpdate = updateResidual / oldWNorm;
            elseif updateResidual == 0
                relativeUpdate = 0;
            else
                relativeUpdate = Inf;
            end
            % Accept the candidate: f^{N} now is saved
            W = Wnew;
            nUpdates = nUpdates + 1;
            % Applying convergence criterion
            meetsTolerance = updateResidual <= cfg.relTol * oldWNorm;
            if meetsTolerance
                normalizedResidual = relativeUpdate;
                status = 'converged';
                break;
            end

            % ----- Update progress recorder -----
            elapsedSeconds = toc(progressTimer);
            if updateIndex == 1 || elapsedSeconds - lastProgressTime >= progressInterval ...
                    || updateIndex == cfg.maxIter
                if hasProgressQueue
                    stopRatio = criterion_ratio(updateResidual, ...
                        cfg.relTol * oldWNorm);
                    report_progress(cfg, 'iteration', model, nUpdates, stopRatio, NaN);
                else
                    iterationsPerSecond = updateIndex / elapsedSeconds;
                    etaSeconds = (cfg.maxIter - updateIndex) / iterationsPerSecond;
                    fprintf(['iter %d/%d | relChange %.3e (tol %.1e) | ' ...
                         'elapsed %.1fs | %.1f iter/s | ETA(maxIter) %.1fmin\n'], ...
                    nUpdates, cfg.maxIter, relativeUpdate, cfg.relTol, ...
                    elapsedSeconds, iterationsPerSecond, etaSeconds / 60);
                end
                lastProgressTime = elapsedSeconds;
            end
        end

        % ----- Unverified exit paths do not yet provide a final residual.
        if ~strcmp(status, 'converged')
            normalizedResidual = NaN;
        end
        result = pack_result(W, status, nUpdates, normalizedResidual, 'sampled_fixed_point');
        result.lastUpdateNorm = updateResidual;
        result.lastRelativeUpdate = relativeUpdate;
        result.iterationSeconds = toc(progressTimer);
        if ~hasProgressQueue
            fprintf(['%s | updates %d | checkedResidual %.3e | ' ...
                     'iteration time %.3fs\n'], ...
                result.status, result.nUpdates, ...
                result.normalizedResidual, result.iterationSeconds);
        end
    else % tau == 0
        integrationTimer = tic;
        lastProgressTime = 0;
        % ----- Check if the initial value already stable
        [dWinitial, failure, initialScale] = compute_dW(W, model);
        if ~isempty(failure)
            error(['riccati:' failure], '%s', failure);
        end

        initialResidual = norm(dWinitial, 'fro') / initialScale;

        if initialResidual <= cfg.relTol
            result = pack_result(W, 'converged', 0, ...
                initialResidual, 'continuous_riccati');
            result.finalTime = 0;
            result.iterationSeconds = toc(integrationTimer);
            return;
        end
        % ----- Setting ode's config -----
        odeOptions = odeset( ...
        'RelTol', cfg.odeRelTol, ...
        'AbsTol', cfg.odeAbsTol, ...
        'Refine', 1);
        odeOptions = odeset(odeOptions, 'Events', ...
        @(T, Wvector) continuous_events(T, Wvector, model, cfg));
        if hasProgressQueue
            odeOptions = odeset(odeOptions, 'OutputFcn', @report_ode_progress);
        end

        % Using compute_dW_vector to find dW and stack in a vector shape
        odeFunction = @(T, Wvector) compute_dW_vector(T, Wvector, model);

        % Load the odeSolver from config
        odeSolver = str2func(cfg.odeSolver);
        [Tvalues, Wtrajectory, ~, ~, eventIndex] = odeSolver(odeFunction, [0, cfg.odeMaxTime], W(:), odeOptions);

        Wend = reshape(Wtrajectory(end, :).', n, n);
        Wend = 0.5 * Wend + 0.5 * Wend.';
        [dWend, failure, residualScale] = compute_dW(Wend, model);
        if ~isempty(failure)
            error(['riccati:' failure], '%s', failure);
        end

        normalizedResidual = norm(dWend, 'fro') / residualScale;

        if any(eventIndex == 1) || norm(Wend, 'fro') >= cfg.maxWNorm
            status = 'inconclusive_W_norm_limit';
        elseif normalizedResidual <= cfg.relTol
            status = 'converged';
        elseif Tvalues(end) < cfg.odeMaxTime * (1 - 1e-12)
            status = 'inconclusive_ode_partial';
        else
            status = 'inconclusive_max_time';
        end

        nUpdates = numel(Tvalues) - 1;

        result = pack_result(Wend, status, nUpdates, ...
            normalizedResidual, 'continuous_riccati');
        result.finalTime = Tvalues(end);
        result.iterationSeconds = toc(integrationTimer);

    end

    function stop = report_ode_progress(T, Wvalues, flag)
        stop = 0;
        if ~isempty(flag)
            return;
        end
        elapsed = toc(integrationTimer);
        if elapsed - lastProgressTime < cfg.progressInterval
            return;
        end
        lastProgressTime = elapsed;
        currentW = reshape(Wvalues(:, end), n, n);
        currentW = 0.5 * currentW + 0.5 * currentW.';
        [currentDW, progressFailure, scale] = compute_dW(currentW, model);
        if isempty(progressFailure)
            stopRatio = criterion_ratio(norm(currentDW, 'fro'), 0.5 * cfg.relTol * scale);
            report_progress(cfg, 'ode', model, 0, stopRatio, T(end));
        end
    end
end

function ratio = criterion_ratio(residual, threshold)
    if threshold > 0
        ratio = residual / threshold;
    elseif residual == 0
        ratio = 0;
    else
        ratio = Inf;
    end
end

function report_progress(cfg, event, model, step, stopRatio, odeTime)
    send(cfg.progressQueue, struct('event', event, ...
        'jobIndex', cfg.progressJobIndex, 'sigmaSquared', model.sigma^2, ...
        'step', step, 'stopRatio', stopRatio, 'odeTime', odeTime));
end

function result = pack_result(W, status, nUpdates, ...
                              normalizedResidual, method)
    result.Wfinal = W;
    result.converged = strcmp(status, 'converged');
    result.status = status;
    result.nUpdates = nUpdates;
    result.normalizedResidual = normalizedResidual;
    result.method = method;
    % Discrete update diagnostics are filled by the sampled branch.
    result.lastUpdateNorm = NaN;
    result.lastRelativeUpdate = NaN;
end

function [value, isterminal, direction] = continuous_events( ...
        ~, Wvector, model, cfg)
    n = size(model.A, 1);
    W = reshape(Wvector, n, n);
    W = 0.5 * W + 0.5 * W.';

    [dW, failure, residualScale] = compute_dW(W, model);
    if ~isempty(failure)
        error(['riccati:' failure], '%s', failure);
    end

    convergenceValue = norm(dW, 'fro') - 0.5 * cfg.relTol * residualScale;

    value = [cfg.maxWNorm - norm(W, 'fro'); convergenceValue];

    isterminal = [1; 1];
    direction = [-1; -1];
end

function dWvector = compute_dW_vector(~, Wvector, model)
    n = size(model.A, 1);
    W = reshape(Wvector, n, n);
    W = 0.5 * W + 0.5 * W.';
    % Compute dW
    [dW, failure] = compute_dW(W, model);
    if ~isempty(failure)
        error(['riccati:' failure], '%s', failure);
    end
    % Put dW into a vector shape
    dWvector = dW(:);
end

function [dW, failure, residualScale] = compute_dW(W, model)
    dW = [];
    failure = '';
    residualScale = NaN;

    A = model.A;
    B = model.B;
    L = model.L;
    R = model.R;
    sigmaSquared = model.sigma^2;

    WB = W * B;
    RplusSigmaSquaredBWB = R + sigmaSquared * (B.' * WB);
    RplusSigmaSquaredBWB = 0.5 * RplusSigmaSquaredBWB + 0.5 * RplusSigmaSquaredBWB.';

    if any(~isfinite(RplusSigmaSquaredBWB(:)))
        failure = 'inconclusive_RplusSigmaSquaredBWB_not_finite';
        return;
    end
    [~, cholFlag] = chol(RplusSigmaSquaredBWB);
    if cholFlag ~= 0
        failure = 'inconclusive_RplusSigmaSquaredBWB_not_positive_definite';
        return;
    end

    controlTerm = WB * (RplusSigmaSquaredBWB \ WB.');
    ATW = A.' * W;
    WA = W * A;
    dW = ATW + WA - controlTerm + L;
    dW = 0.5 * dW + 0.5 * dW.';

    if any(~isfinite(dW(:)))
        dW = [];
        failure = 'inconclusive_dW_not_finite';
        return;
    end
    residualScale = max(1, ...
        norm(ATW, 'fro') + norm(WA, 'fro') + ...
        norm(controlTerm, 'fro') + norm(L, 'fro'));

    if ~isfinite(residualScale)
        dW = [];
        residualScale = NaN;
        failure = 'inconclusive_residual_scale_not_finite';
        return;
    end
end

function [Wnew, failure] = compute_f(W, model, precomputed, maxWNorm)
    Wnew = [];
    failure = '';
    n = size(W, 1);

    % Loading variables from model and precomputed set
    R = model.R;
    sigmaSquared = model.sigma^2;
    tau = model.tau;
    L = model.L;
    expAsBGrid = precomputed.expAsBGrid;
    quadratureWeights = precomputed.quadratureWeights;
    expATau = precomputed.expATau;

    % Compute M and check if there's any failure before proceding
    [M, Mfailure] = compute_M(W, expAsBGrid, R, sigmaSquared, quadratureWeights);
    if ~isempty(Mfailure)
        failure = Mfailure;
        return;
    end

    % Compute W^(1/2) and check if there's failure in computing sqrtW before proceeding
    [sqrtW, ~] = sqrtm(W);
    if ~isreal(sqrtW)
        failure = 'inconclusive_sqrtW_not_real';
        return;
    end
    if any(~isfinite(sqrtW(:)))
        failure = 'inconclusive_sqrtW_not_finite';
        return
    end

    % compute I + W^(1/2)MW^(1/2) and check if IplusSqrtWMSqrtW is finite and positive definite
    IplusSqrtWMSqrtW = eye(n) + sqrtW * M * sqrtW;
    IplusSqrtWMSqrtW = 0.5 * IplusSqrtWMSqrtW + 0.5 * IplusSqrtWMSqrtW.';
    if any(~isfinite(IplusSqrtWMSqrtW(:)))
        failure = 'inconclusive_IplusSqrtWMSqrtW_not_finite';
        return;
    end
    [~, cholFlag] = chol(IplusSqrtWMSqrtW);
    if cholFlag ~= 0
        failure = 'inconclusive_IplusSqrtWMSqrtW_not_positive_definite';
        return;
    end

    % Compute f(W) and validate the candidate before accepting it.
    Wnew = expATau.' * sqrtW * (IplusSqrtWMSqrtW \ sqrtW) * expATau + tau * L;
    if any(~isfinite(Wnew(:)))
        failure = 'inconclusive_W_not_finite';
        return;
    end
    Wnew = 0.5 * Wnew + 0.5 * Wnew.';
    if norm(Wnew, 'fro') >= maxWNorm
        failure = 'inconclusive_W_norm_limit';
        return;
    end

end

function [M, failure] = compute_M(W, expAsBGrid, R, sigmaSquared, quadratureWeights)
    n = size(W, 1);
    M = zeros(n);
    failure = '';

    for j = 1:numel(quadratureWeights)
        expAsB = expAsBGrid(:, :, j);

        % Compute Q(s)
        Q = R + sigmaSquared * (expAsB.' * W * expAsB);
        Q = (Q + Q.') / 2; % prevent loss of symmetry due to numerical error

        % Check if Q is finite and positive definite
        if any(~isfinite(Q(:)))
            M = [];
            failure = 'inconclusive_Q_not_finite';
            return;
        end
        [~, cholFlag] = chol(Q);
        if cholFlag ~= 0
            M = [];
            failure = 'inconclusive_Q_not_positive_definite';
            return;
        end

        % Compute the integrand for M
        integrand = expAsB * (Q \ expAsB.');

        % Add to the integral
        M = M + quadratureWeights(j) * integrand;
    end

    M = (M + M.') / 2;
end

function validate_model(model, maxWNorm)
    A = model.A;
    B = model.B;
    L = model.L;
    R = model.R;
    P = model.P;

    n = size(A, 1);
    m = size(B, 2);

    % Check if they are real matrices with neither Nan nor Inf.
    matrixNames = {'A', 'B', 'L', 'R', 'P'};

    for k = 1:numel(matrixNames)
        matrixName = matrixNames{k};
        matrixValue = model.(matrixName);

        if ~isnumeric(matrixValue) || ~isreal(matrixValue) || ...
                ~all(isfinite(matrixValue(:)))

            error('riccati:InvalidModel', ...
                'model.%s must contain finite real numeric values.', matrixName);
        end
    end

    % Dimension check
    if ~isequal(size(A), [n, n])
        error('riccati:InvalidModel', ...
            'model.A must be n-by-n.');
    end
    if m < 1 || ~isequal(size(B), [n, m])
        error('riccati:InvalidModel', ...
            'model.B must be n-by-m with m >= 1.');
    end
    if ~isequal(size(R), [m, m])
        error('riccati:InvalidModel', ...
            'model.R must be m-by-m.');
    end
    if ~isequal(size(L), [n, n])
        error('riccati:InvalidModel', ...
            'model.L must be n-by-n.');
    end
    if ~isequal(size(P), [n, n])
        error('riccati:InvalidModel', ...
            'model.P must be n-by-n.');
    end

    % Check if L, R, P are symmetric.
    symmetricMatrixNames = {'L', 'R', 'P'};
    symmetryTolerance = 1e-12;

    for k = 1:numel(symmetricMatrixNames)
        matrixName = symmetricMatrixNames{k};
        matrixValue = model.(matrixName);

        symmetryError = norm(matrixValue - matrixValue.', 'fro');
        matrixScale = max(1, norm(matrixValue, 'fro'));

        if symmetryError > symmetryTolerance * matrixScale
            error('riccati:InvalidModel', ...
                'model.%s must be symmetric.', matrixName);
        end
    end

    % Check if L, R are positive definite
    positiveMatrixNames = {'L', 'R'};

    for k = 1:numel(positiveMatrixNames)
        matrixName = positiveMatrixNames{k};
        matrixValue = model.(matrixName);
        matrixValue = (matrixValue + matrixValue.') / 2;

        [~, cholFlag] = chol(matrixValue);

        if cholFlag ~= 0
            error('riccati:InvalidModel', ...
                'model.%s must be positive definite.', matrixName);
        end
    end

    % Check if P is positive semidefinite
    symmetricP = (P + P.') / 2;
    psdTolerance = 1e-12 * max(1, norm(symmetricP, 2));

    if min(eig(symmetricP)) < -psdTolerance
        error('riccati:InvalidModel', ...
            'model.P must be positive semidefinite.');
    end


    % Check if P is below cfg.maxWNorm
    if norm(P, 'fro') >= maxWNorm
        error('riccati:InvalidInitialMatrix', ...
            'The norm of model.P must be below cfg.maxWNorm.');
    end

    % Check if tau and sigma are valid finite real scalars both are >= 0
    scalarNames = {'tau', 'sigma'};

    for k = 1:numel(scalarNames)
        scalarName = scalarNames{k};
        scalarValue = model.(scalarName);

        if ~isnumeric(scalarValue) || ~isreal(scalarValue) || ...
                ~isscalar(scalarValue) || ~isfinite(scalarValue) || ...
                scalarValue < 0

            error('riccati:InvalidModel', ...
                'model.%s must be a finite real scalar >= 0.', scalarName);
        end
    end


end

function validate_config(cfg, tau)
    requiredFields = {'relTol', 'maxWNorm'};
    if tau > 0
        requiredFields = [requiredFields, {'Ns', 'maxIter'}];
    else
        requiredFields = [requiredFields, ...
            {'odeSolver', 'odeMaxTime', 'odeRelTol', 'odeAbsTol'}];
    end
    for k = 1:numel(requiredFields)
        fieldName = requiredFields{k};
        % Check if num of fields are correct
        if ~isfield(cfg, fieldName)
            error('riccati:MissingConfig', ...
                'Missing cfg.%s.', fieldName);
        end
    end

    % Check if relTol is valid input
    if ~isnumeric(cfg.relTol) || ~isreal(cfg.relTol) || ...
    ~isscalar(cfg.relTol) || ~isfinite(cfg.relTol) || ...
    cfg.relTol <= 0
        error('riccati:InvalidConfig', ...
            'cfg.relTol must be a finite real number > 0.');
    end
    % Check if maxWNorm is valid input
    if ~isnumeric(cfg.maxWNorm) || ~isreal(cfg.maxWNorm) || ...
    ~isscalar(cfg.maxWNorm) || ~isfinite(cfg.maxWNorm) || ...
    cfg.maxWNorm <= 0
        error('riccati:InvalidConfig', ...
            'cfg.maxWNorm must be a finite real number > 0.');
    end

    if tau > 0
         % Check if Ns is valid input
        if ~isnumeric(cfg.Ns) || ~isreal(cfg.Ns) || ...
            ~isscalar(cfg.Ns) || ~isfinite(cfg.Ns) || cfg.Ns < 2 ||...
            cfg.Ns ~= fix(cfg.Ns)
            error('riccati:InvalidConfig', ...
                'cfg.Ns must be a finite integer >=  2.');
        end
        % Check if maxIter is valid input
        if ~isnumeric(cfg.maxIter) || ~isreal(cfg.maxIter) || ...
        ~isscalar(cfg.maxIter) || ~isfinite(cfg.maxIter) || ...
        cfg.maxIter < 1 || cfg.maxIter ~= fix(cfg.maxIter)
            error('riccati:InvalidConfig', ...
                'cfg.maxIter must be a finite real integer >= 1.');
        end
    else
        % Check ODE parameters.
        positiveFields = {'odeMaxTime', 'odeRelTol', 'odeAbsTol'};

        for k = 1:numel(positiveFields)
            fieldName = positiveFields{k};
            value = cfg.(fieldName);
            % Check if the entry is valid input
            if ~isnumeric(value) || ~isreal(value) || ...
                    ~isscalar(value) || ~isfinite(value) || value <= 0
                error('riccati:InvalidConfig', ...
                    'cfg.%s must be a finite real scalar > 0.', fieldName);
            end
        end

        % Check if odeSolver is the ones we accepted
        solverName = cfg.odeSolver;
        isSingleName = (ischar(solverName) && isrow(solverName)) || ...
                       (isstring(solverName) && isscalar(solverName));
        if ~isSingleName || ...
                ~any(strcmp(solverName, {'ode15s', 'ode45'}))
            error('riccati:InvalidConfig', ...
                'cfg.odeSolver must be ode15s or ode45.');

        end
    end

end
