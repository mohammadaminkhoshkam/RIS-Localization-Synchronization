function save_crlb_outputs(outputFolder, results, summaryTable, ...
    parameterTable, priorTable, figures)
% Save tables, figures, and MATLAB results.

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

writetable(summaryTable, fullfile(outputFolder, 'summary_bounds.csv'));
writetable(parameterTable, ...
    fullfile(outputFolder, 'all_parameter_bounds.csv'));
writetable(results.sweepTable, ...
    fullfile(outputFolder, 'snr_sweep_bounds.csv'));

if ~isempty(priorTable)
    writetable(priorTable, ...
        fullfile(outputFolder, 'prior_sensitivity.csv'));
end

figureFields = {'geometry', 'mainBounds', 'ueComponents', ...
    'risComponents', 'synchronization', 'defaultBars', ...
    'fimDiagnostics'};
figureNames = {'01_geometry.png', '02_main_bounds.png', ...
    '03_ue_components.png', '04_ris_components.png', ...
    '05_sync_beta.png', '06_default_bars.png', ...
    '07_fim_diagnostics.png'};

if isfield(figures, 'priorSensitivity')
    figureFields{end+1} = 'priorSensitivity';
    figureNames{end+1} = '08_prior_sensitivity.png';
end

for i = 1:numel(figureFields)
    fig = figures.(figureFields{i});
    if isgraphics(fig)
        save_project_figure(fig, fullfile(outputFolder, figureNames{i}));
    end
end

resultsForFile = results;
save(fullfile(outputFolder, 'crlb_results.mat'), ...
    'resultsForFile', 'summaryTable', 'priorTable', '-v7');
fprintf('Results saved in: %s\n', outputFolder);
end


function save_project_figure(fig, filename)
% Save one figure with a fallback method.
try
    exportgraphics(fig, filename, 'Resolution', 180);
catch
    saveas(fig, filename);
end
end
