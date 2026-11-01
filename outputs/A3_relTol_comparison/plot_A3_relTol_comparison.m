%% Compare three relative tolerances for A_3 and the two input matrices
% This script reads the adjacent results.mat and recreates both figures.
% sigmaCurveI and sigmaCurveOnes contain sigma, with one column per relTol.

comparisonDir = fileparts(mfilename('fullpath'));
data = load(fullfile(comparisonDir, 'results.mat'));

[figI, axI] = plot_comparison(data.tauSet{1}, data.sigmaCurveI, data.relTols, 2);
[figOnes, axOnes] = plot_comparison(data.tauSet{2}, data.sigmaCurveOnes, data.relTols, 1);

savefig(figI, fullfile(comparisonDir, 'A3_B_I5_relTol_comparison.fig'));
exportgraphics(axI, fullfile(comparisonDir, 'A3_B_I5_relTol_comparison.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
savefig(figOnes, fullfile(comparisonDir, 'A3_B_ones5_relTol_comparison.fig'));
exportgraphics(axOnes, fullfile(comparisonDir, 'A3_B_ones5_relTol_comparison.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');

fprintf('A_3 comparison figures saved to:\n%s\n', comparisonDir);

function [fig, ax] = plot_comparison(tauValues, sigmaCurves, relTols, xTickStep)
    tauValues = tauValues(:);
    assert(size(sigmaCurves, 1) == numel(tauValues));
    assert(size(sigmaCurves, 2) == numel(relTols));

    sigmaSquared = sigmaCurves.^2;
    sigmaSquared(~isfinite(sigmaSquared) | sigmaSquared <= 0) = NaN;
    positiveValues = sigmaSquared(isfinite(sigmaSquared));
    assert(~isempty(positiveValues), 'No positive boundary values to plot.');

    colors = [0.0000 0.4470 0.7410; ...
              0.8500 0.3250 0.0980; ...
              0.4940 0.1840 0.5560];
    lineStyles = {'-', '--', ':'};
    markers = {'o', 's', '^'};

    fig = figure('Color', 'w', 'Units', 'pixels', ...
        'Position', [100 100 968 660]);
    ax = axes('Parent', fig, 'Position', [0.135 0.150 0.835 0.810]);
    handles = semilogy(ax, tauValues, sigmaSquared, 'LineWidth', 1.7);
    for k = 1:numel(relTols)
        % Stagger marker locations to distinguish coincident curves.
        set(handles(k), 'Color', colors(k, :), ...
            'LineStyle', lineStyles{k}, 'Marker', markers{k}, ...
            'MarkerIndices', k:3:numel(tauValues), 'MarkerSize', 6.5, ...
            'MarkerFaceColor', colors(k, :), 'MarkerEdgeColor', colors(k, :), ...
            'DisplayName', sprintf('relTol = 10^{%d}', round(log10(relTols(k)))));
    end

    firstExponent = floor(log10(min(positiveValues)));
    lastExponent = ceil(log10(max(positiveValues)));
    if firstExponent == lastExponent
        lastExponent = firstExponent + 1;
    end
    set(ax, 'YScale', 'log', ...
        'XLim', [0 max(tauValues)], 'XTick', 0:xTickStep:max(tauValues), ...
        'YLim', 10.^[firstExponent lastExponent], ...
        'YTick', 10.^(firstExponent:lastExponent), ...
        'FontName', 'Helvetica', 'FontSize', 20, ...
        'LineWidth', 1.1, 'TickDir', 'in', 'Box', 'on', ...
        'XGrid', 'on', 'YGrid', 'on', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'on', ...
        'GridColor', [0 0 0], 'MinorGridColor', [0 0 0], ...
        'GridAlpha', 0.24, 'MinorGridAlpha', 0.16, ...
        'GridLineStyle', ':', 'MinorGridLineStyle', ':', 'Layer', 'top');
    xlabel(ax, '\tau', 'Interpreter', 'tex', 'FontSize', 22);
    ylabel(ax, '\sigma^2_{max}', 'Interpreter', 'tex', 'FontSize', 22);
    legendHandle = legend(ax, handles, 'Location', 'northeast', ...
        'Interpreter', 'tex', 'FontSize', 18, 'LineWidth', 0.8, ...
        'Box', 'on', 'Color', 'w', 'EdgeColor', 'k');
    legendHandle.IconColumnWidth = 42;
end
