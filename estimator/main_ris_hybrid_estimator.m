%% Hybrid MAP/ML estimator with variable projection

clearvars;
clear functions;
close all;
clc;
format compact;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptFolder = pwd;
else
    scriptFolder = fileparts(scriptPath);
end
addpath(scriptFolder);

%% 1. Configuration

cfg = ris_estimator_config();
K = size(cfg.true.p0, 2);
idx = ris_parameter_indices(K);
state = ris_build_states(cfg, idx);

configuration = table(cfg.N, cfg.L, cfg.deltaF, cfg.snapshotInterval, ...
    cfg.snrDb, K, idx.nPhysical, cfg.optimizer.numStarts, ...
    'VariableNames', {'N', 'L', 'DeltaF_Hz', 'SnapshotInterval_s', ...
    'SNR_dB', 'K', 'PhysicalDimension', 'NumStarts'})

%% 2. Simulate the received OFDM samples

data = ris_simulate_data(state.true, cfg.true.beta, cfg, idx);

signalSummary = table(data.signalPower, data.sigma2, ...
    10*log10(data.signalPower/data.sigma2), ...
    'VariableNames', {'SignalPower', 'NoiseVariance', 'MeasuredSNR_dB'})

%% 3. Run the concentrated hybrid MAP/ML estimator

tic;
result = ris_estimate_hybrid_map(data, state, cfg, idx);
elapsedTime = toc;

startSummary = table((1:cfg.optimizer.numStarts).', result.startCosts, ...
    result.finalCosts, result.exitFlags, ...
    'VariableNames', {'Start', 'InitialObjective', ...
    'FinalObjective', 'ExitFlag'})

syncInitialization = table( ...
    result.syncInitialization.initialCost, ...
    result.syncInitialization.bestCost, ...
    result.syncInitialization.initialEtaPpm, ...
    result.syncInitialization.bestEtaPpm, ...
    result.syncInitialization.initialDeltaNs, ...
    result.syncInitialization.bestDeltaNs, ...
    'VariableNames', {'InitialObjective', 'GridObjective', ...
    'InitialEta_ppm', 'GridEta_ppm', ...
    'InitialDelta_ns', 'GridDelta_ns'})

fprintf('Best start: %d\n', result.bestStart);
fprintf('Elapsed time: %.3f seconds\n', elapsedTime);
fprintf('Used fmincon: %d\n', result.usedFmincon);
normalizedDataCost = result.dataCost/numel(data.y);
fprintf('Normalized data objective: %.4f\n', normalizedDataCost);
if normalizedDataCost > cfg.optimizer.maxNormalizedDataCost
    warning(['The selected solution has a large normalized data cost. ' ...
        'Increase the synchronization search span or the number of starts.']);
end

%% 4. Local information and identifiability

info = ris_information_diagnostics(result, data, state, cfg, idx);

fimSummary = table(info.rank, info.dimension, info.conditionNumber, ...
    info.fullRank, 'VariableNames', {'Rank', 'Dimension', ...
    'ConditionNumber', 'FullRank'})

%% 5. Estimation tables

tables = ris_result_tables(result, state, info, cfg, idx);

summaryResults = tables.summary
physicalParameterResults = tables.physical
channelResults = tables.channel

%% 6. Figures

figures = ris_plot_results(result, state, data, info, cfg, idx);

%% 7. Save outputs

if cfg.saveResults
    outputFolder = fullfile(scriptFolder, cfg.outputFolder);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    writetable(tables.summary, fullfile(outputFolder, 'summary.csv'));
    writetable(tables.physical, ...
        fullfile(outputFolder, 'physical_parameters.csv'));
    writetable(tables.channel, ...
        fullfile(outputFolder, 'channel_parameters.csv'));
    writetable(startSummary, ...
        fullfile(outputFolder, 'multi_start_summary.csv'));
    writetable(syncInitialization, ...
        fullfile(outputFolder, 'sync_initialization.csv'));

    resultToSave = result;
    save(fullfile(outputFolder, 'hybrid_map_estimate.mat'), ...
        'cfg', 'state', 'data', 'resultToSave', 'info', '-v7');

    try
        exportgraphics(figures.geometry, ...
            fullfile(outputFolder, 'geometry.png'), 'Resolution', 180);
        exportgraphics(figures.diagnostics, ...
            fullfile(outputFolder, 'diagnostics.png'), 'Resolution', 180);
        if isfield(figures, 'initialization')
            exportgraphics(figures.initialization, ...
                fullfile(outputFolder, 'sync_initialization.png'), ...
                'Resolution', 180);
        end
    catch
        saveas(figures.geometry, fullfile(outputFolder, 'geometry.png'));
        saveas(figures.diagnostics, fullfile(outputFolder, 'diagnostics.png'));
        if isfield(figures, 'initialization')
            saveas(figures.initialization, ...
                fullfile(outputFolder, 'sync_initialization.png'));
        end
    end
    fprintf('Results saved in: %s\n', outputFolder);
end
