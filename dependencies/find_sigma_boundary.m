function boundary = find_sigma_boundary(model, solverCfg, searchCfg)

    sigmaSquareLow = 0;
    modelSearch = model;
    modelSearch.sigma = sqrt(sigmaSquareLow);
    lowProbe = riccati(modelSearch, solverCfg);

    sigmaSquareHigh = 4;
    modelSearch.sigma = sqrt(sigmaSquareHigh);
    highProbe = riccati(modelSearch, solverCfg);

    hasBracket = lowProbe.converged && ~highProbe.converged;

    if hasBracket

        % Bisect the noise variance sigma^2 for one fixed (A, B, tau).
        while sigmaSquareHigh - sigmaSquareLow > ...
        searchCfg.sigmaSquaredRelTol * max(sigmaSquareHigh, 1e-12)
            sigmaSquareMid = (sigmaSquareLow + sigmaSquareHigh) / 2;
            modelSearch.sigma = sqrt(sigmaSquareMid);

            midProbe = riccati(modelSearch, solverCfg);

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
end
