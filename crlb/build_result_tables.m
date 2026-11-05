function [statusTable, summaryTable, parameterTable] = ...
    build_result_tables(results)
% Build the main result tables.

R = results.default;
K = size(results.config.risPositionM, 2);

statusTable = table( ...
    R.snrDb, R.noiseVariance, R.infoEquivalent.rank, ...
    R.infoEquivalent.dimension, R.infoEquivalent.conditionNumber, ...
    'VariableNames', {'SNR_dB', 'NoiseVariance', 'FIM_Rank', ...
    'FIM_Dimension', 'ScaledFIM_ConditionNumber'});

metric = {'UE position error bound'; 'UE velocity error bound'};
value = [R.uePEB_m; R.ueVEB_mps];
unit = {'m'; 'm/s'};

for k = 1:K
    metric{end+1, 1} = sprintf('RIS %d position error bound', k);
    value(end+1, 1) = R.risPEB_m(k);
    unit{end+1, 1} = 'm';

    metric{end+1, 1} = sprintf('RIS %d velocity error bound', k);
    value(end+1, 1) = R.risVEB_mps(k);
    unit{end+1, 1} = 'm/s';
end

metric{end+1, 1} = 'Relative frequency error eta';
value(end+1, 1) = R.etaStdPpm;
unit{end+1, 1} = 'ppm';

metric{end+1, 1} = 'Clock offset delta_c';
value(end+1, 1) = R.deltaStdNs;
unit{end+1, 1} = 'ns';

summaryTable = table(metric, value, unit, ...
    'VariableNames', {'EstimatedQuantity', 'StdLowerBound', 'Unit'});
parameterTable = results.parameterTable;
end
