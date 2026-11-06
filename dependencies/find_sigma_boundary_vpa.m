function boundary = find_sigma_boundary_vpa(model, solverCfg, searchCfg)

    probeIndex = 0;
    nBisections = 0;
    showProgress = ~isfield(solverCfg, 'progressQueue');

    sigmaSquareLow = 0;
    modelSearch = model;
    modelSearch.sigma = sqrt(sigmaSquareLow);
    probeIndex = probeIndex + 1;
    lowProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

    sigmaSquareHigh = 4;
    modelSearch.sigma = sqrt(sigmaSquareHigh);
    probeIndex = probeIndex + 1;
    highProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

    hasBracket = lowProbe.converged && ~highProbe.converged;

    if hasBracket

        % Bisect the noise variance sigma^2 for one fixed (A, B, tau).
        while sigmaSquareHigh - sigmaSquareLow > ...
        searchCfg.sigmaSquaredRelTol * max(sigmaSquareHigh, 1e-12)
            nBisections = nBisections + 1;
            if showProgress
                relativeWidth = (sigmaSquareHigh - sigmaSquareLow) / max(sigmaSquareHigh, 1e-12);
                fprintf(['Bisection %d | sigma^2 bracket [%.6g, %.6g]' ...
                    ' | relative width %.2f%% (target %.2f%%)\n'], ...
                    nBisections, sigmaSquareLow, sigmaSquareHigh, ...
                    100 * relativeWidth, 100 * searchCfg.sigmaSquaredRelTol);
            end
            sigmaSquareMid = (sigmaSquareLow + sigmaSquareHigh) / 2;
            modelSearch.sigma = sqrt(sigmaSquareMid);

            probeIndex = probeIndex + 1;
            midProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

            if midProbe.converged
                sigmaSquareLow = sigmaSquareMid;
                lowProbe = midProbe;
            else
                sigmaSquareHigh = sigmaSquareMid;
                highProbe = midProbe;
            end
        end


    else
        sigmaSquareLow = NaN;
        sigmaSquareHigh = NaN;
    end
    boundary.tau = modelSearch.tau;
    boundary.sigmaSquareLow = sigmaSquareLow;
    boundary.sigmaSquareHigh = sigmaSquareHigh;
    boundary.lowProbe = lowProbe;
    boundary.highProbe = highProbe;
    if showProgress
        if hasBracket
            relativeWidth = (sigmaSquareHigh - sigmaSquareLow) / max(sigmaSquareHigh, 1e-12);
            fprintf(['Boundary search finished: %d probes, %d bisections' ...
                ' | relative width %.2f%%\n'], ...
                probeIndex, nBisections, 100 * relativeWidth);
        else
            fprintf('Boundary search stopped after %d probes: invalid initial bracket.\n', probeIndex);
        end
    end
end

function probe = run_sigma_probe(model, solverCfg, probeIndex)
    if ~isfield(solverCfg, 'progressQueue')
        fprintf('\nProbe %d | sigma = %.8g | sigma^2 = %.8g\n', ...
            probeIndex, model.sigma, model.sigma^2);
    end
    probe = riccati_vpa(model, solverCfg);
end
