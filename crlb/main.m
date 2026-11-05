%% CRLB simulation
clearvars;
clear functions;
close all;
clc;
format compact;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    projectFolder = pwd;
else
    projectFolder = fileparts(scriptPath);
end
addpath(projectFolder);


cfg = create_crlb_config();
runPriorSensitivity = true;
saveOutputs = true;

%% Compute bounds
results = simulate_all_crlb_bounds(cfg);
[calculationStatus, summaryTable, allParameterBounds] = ...
    build_result_tables(results);

calculationStatus
summaryTable
allParameterBounds

%% Plot results
figures = plot_crlb_bounds(results);
figures.geometry = plot_system_geometry(results.config);
figures.fimDiagnostics = plot_fim_diagnostics(results);

%% Analyze RIS prior
if runPriorSensitivity
    [priorSensitivityTable, figures.priorSensitivity] = ...
        analyze_ris_prior(results.config);
    priorSensitivityTable
else
    priorSensitivityTable = table();
end

print_crlb_summary(results.default);

%% Save outputs
if saveOutputs
    outputFolder = fullfile(projectFolder, 'crlb_results');
    save_crlb_outputs(outputFolder, results, summaryTable, ...
        allParameterBounds, priorSensitivityTable, figures);
end
