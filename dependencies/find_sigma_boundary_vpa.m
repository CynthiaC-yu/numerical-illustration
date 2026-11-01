function boundary = find_sigma_boundary_vpa(model, solverCfg, searchCfg)

    probeIndex = 0;
    nBisections = 0;
    showProgress = ~isfield(solverCfg, 'progressQueue');

    sigmaLow = 0;
    modelSearch = model;
    modelSearch.sigma = sigmaLow;
    probeIndex = probeIndex + 1;
    lowProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

    sigmaHigh = 2;
    modelSearch.sigma = sigmaHigh;
    probeIndex = probeIndex + 1;
    highProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

    hasBracket = lowProbe.converged && ~highProbe.converged;

    if hasBracket

        % Binary search loop for one fixed (A, B, tau)
        while sigmaHigh^2 - sigmaLow^2 > ...
        searchCfg.sigmaSquaredRelTol * max(sigmaHigh^2, 1e-12)
            nBisections = nBisections + 1;
            if showProgress
                relativeWidth = (sigmaHigh^2 - sigmaLow^2) / max(sigmaHigh^2, 1e-12);
                fprintf(['Bisection %d | sigma^2 bracket [%.6g, %.6g]' ...
                    ' | relative width %.2f%% (target %.2f%%)\n'], ...
                    nBisections, sigmaLow^2, sigmaHigh^2, ...
                    100 * relativeWidth, 100 * searchCfg.sigmaSquaredRelTol);
            end
            sigmaMid = (sigmaLow + sigmaHigh) / 2;
            modelSearch.sigma = sigmaMid;

            probeIndex = probeIndex + 1;
            midProbe = run_sigma_probe(modelSearch, solverCfg, probeIndex);

            if midProbe.converged
                sigmaLow = sigmaMid;
                lowProbe = midProbe;
            else
                sigmaHigh = sigmaMid;
                highProbe = midProbe;
            end
        end


    else
        sigmaLow = NaN;
        sigmaHigh = NaN;
    end
    boundary.tau = modelSearch.tau;
    boundary.sigmaLow = sigmaLow;
    boundary.sigmaHigh = sigmaHigh;
    boundary.lowProbe = lowProbe;
    boundary.highProbe = highProbe;
    if showProgress
        if hasBracket
            relativeWidth = (sigmaHigh^2 - sigmaLow^2) / max(sigmaHigh^2, 1e-12);
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
