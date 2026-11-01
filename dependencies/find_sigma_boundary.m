function boundary = find_sigma_boundary(model, solverCfg, searchCfg)

    sigmaLow = 0;
    modelSearch = model;
    modelSearch.sigma = sigmaLow;
    lowProbe = riccati(modelSearch, solverCfg);

    sigmaHigh = 2;
    modelSearch.sigma = sigmaHigh;
    highProbe = riccati(modelSearch, solverCfg);

    hasBracket = lowProbe.converged && ~highProbe.converged;

    if hasBracket

        % Binary search loop for one fixed (A, B, tau)
        while sigmaHigh^2 - sigmaLow^2 > ...
        searchCfg.sigmaSquaredRelTol * max(sigmaHigh^2, 1e-12)
            sigmaMid = (sigmaLow + sigmaHigh) / 2;
            modelSearch.sigma = sigmaMid;

            midProbe = riccati(modelSearch, solverCfg);

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
end
