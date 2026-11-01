function result = riccati_vpa(model, cfg)
%RICCATI_VPA Solve the sampled or continuous Riccati problem in VPA arithmetic.
% For tau=0, use odeSolver='rk45_vpa' and supply odeMaxStep and ODE tolerances.
% ----- Check if riccati solver received a complete signature.
    if nargin < 2
        error('riccati:MissingInput', ...
            'Both model and cfg must be provided.')
    end

    % ----- Check if model and config input is complete and valid.
    validate_model(model, cfg.maxWNorm);
    validate_config(cfg, model.tau);

    % Set precision for this solver call.
    previousDigits = digits(cfg.precisionDigits);

    % Restore the previous precision when this call finishes.
    precisionCleanup = onCleanup(@() digits(previousDigits));


    % Convert model parameters to high-precision numbers.
    modelFields = {'A', 'B', 'L', 'R', 'P', 'tau', 'sigma'};

    for fieldIndex = 1:numel(modelFields)
        fieldName = modelFields{fieldIndex};

        model.(fieldName) = vpa( ...
            sym(model.(fieldName), 'f'), cfg.precisionDigits);
    end
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
        h = tau / (Ns - 1);
        sGrid = h * sym(0:Ns-1);

        quadratureWeights = h * sym(ones(1, Ns));
        quadratureWeights(1) = h / 2;
        quadratureWeights(end) = h / 2;

        expAsBGrid = sym(zeros(n, m, Ns));
        for j = 1:Ns
            expAsBGrid(:, :, j) = vpa(expm(A * sGrid(j)) * B);
        end

        expATau = vpa(expm(A * tau));
        % Assigning vars
        precomputed.expAsBGrid = expAsBGrid;
        precomputed.quadratureWeights = quadratureWeights;
        precomputed.expATau = expATau;

        % ----- The loop for updating W -----
        % Initiate key variables
        status = 'inconclusive_max_iter';
        nUpdates = 0;
        normalizedResidual = NaN;
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
            % Normalize with the same scale as the min criterion.
            normalizationScale = min(1, oldWNorm);
            normalizedUpdate = criterion_ratio(updateResidual, normalizationScale);
            % Accept the candidate: f^{N} now is saved
            W = Wnew;
            nUpdates = nUpdates + 1;
            % Applying convergence criterion
            meetsTolerance = updateResidual <= cfg.relTol * normalizationScale;
            if meetsTolerance
                normalizedResidual = normalizedUpdate;
                status = 'converged';
                break;
            end

            % ----- Update progress recorder -----
            elapsedSeconds = toc(progressTimer);
            if updateIndex == 1 || elapsedSeconds - lastProgressTime >= progressInterval ...
                    || updateIndex == cfg.maxIter
                if hasProgressQueue
                    stopRatio = criterion_ratio(updateResidual, ...
                        cfg.relTol * normalizationScale);
                    report_progress(cfg, 'iteration', model, nUpdates, stopRatio, NaN);
                else
                    iterationsPerSecond = updateIndex / elapsedSeconds;
                    etaSeconds = (cfg.maxIter - updateIndex) / iterationsPerSecond;
                    fprintf(['iter %d/%d | normalizedChange %.3e (tol %.1e) | ' ...
                         'elapsed %.1fs | %.1f iter/s | ETA(maxIter) %.1fmin\n'], ...
                    nUpdates, cfg.maxIter, normalizedUpdate, cfg.relTol, ...
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
        result.iterationSeconds = toc(progressTimer);
        if ~hasProgressQueue
            fprintf(['%s | updates %d | checkedResidual %.3e | ' ...
                     'iteration time %.3fs\n'], ...
                result.status, result.nUpdates, ...
                result.normalizedResidual, result.iterationSeconds);
        end
    else % tau == 0
        result = integrate_riccati_vpa(model, cfg, @compute_dW);
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
        'jobIndex', cfg.progressJobIndex, 'sigmaSquared', double(model.sigma^2), ...
        'step', step, 'stopRatio', double(stopRatio), 'odeTime', double(odeTime)));
end

function result = pack_result(W, status, nUpdates, ...
                              normalizedResidual, method)
    result.Wfinal = W;
    result.converged = strcmp(status, 'converged');
    result.status = status;
    result.nUpdates = nUpdates;
    result.normalizedResidual = normalizedResidual;
    result.method = method;
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
    residualScale = min(1, norm(W, 'fro'));

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
    requiredFields = {'relTol', 'maxWNorm', 'precisionDigits', 'maxIter'};
    if tau > 0
        requiredFields = [requiredFields, {'Ns'}];
    else
        requiredFields = [requiredFields, ...
            {'odeSolver', 'odeMaxTime', 'odeMaxStep', 'odeRelTol', 'odeAbsTol'}];
    end
    for k = 1:numel(requiredFields)
        fieldName = requiredFields{k};
        % Check if num of fields are correct
        if ~isfield(cfg, fieldName)
            error('riccati:MissingConfig', ...
                'Missing cfg.%s.', fieldName);
        end
    end

    % Precision and iteration budget apply to both solver branches.
    if ~isnumeric(cfg.precisionDigits) || ~isreal(cfg.precisionDigits) || ...
            ~isscalar(cfg.precisionDigits) || ~isfinite(cfg.precisionDigits) || ...
            cfg.precisionDigits < 16 || cfg.precisionDigits ~= fix(cfg.precisionDigits)
        error('riccati:InvalidConfig', ...
            'cfg.precisionDigits must be an integer >= 16.');
    end
    if ~isnumeric(cfg.maxIter) || ~isreal(cfg.maxIter) || ...
            ~isscalar(cfg.maxIter) || ~isfinite(cfg.maxIter) || ...
            cfg.maxIter < 1 || cfg.maxIter ~= fix(cfg.maxIter)
        error('riccati:InvalidConfig', ...
            'cfg.maxIter must be a finite integer >= 1.');
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
    else
        % Check ODE parameters.
        positiveFields = {'odeMaxTime', 'odeMaxStep', 'odeRelTol', 'odeAbsTol'};

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

        % The continuous branch uses a Runge-Kutta solver in VPA arithmetic.
        solverName = cfg.odeSolver;
        isSingleName = (ischar(solverName) && isrow(solverName)) || ...
                       (isstring(solverName) && isscalar(solverName));
        if ~isSingleName || ...
                ~strcmp(solverName, 'rk45_vpa')
            error('riccati:InvalidConfig', ...
                'For tau=0, cfg.odeSolver must be rk45_vpa.');

        end
    end

end
