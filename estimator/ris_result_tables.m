function tables = ris_result_tables(result, state, info, cfg, idx)
% Create physical, channel, and summary tables.

labels = physical_labels(idx);
units = repmat({'SI'}, idx.nPhysical, 1);
factor = ones(idx.nPhysical, 1);

units(idx.u0) = {'m'; 'm'; 'm'};
units(idx.uDot) = {'m/s'; 'm/s'; 'm/s'};
for k = 1:idx.K
    units(idx.p0(:, k)) = {'m'; 'm'; 'm'};
    units(idx.pDot(:, k)) = {'m/s'; 'm/s'; 'm/s'};
end
units{idx.eta} = 'ppm';
units{idx.deltaC} = 'ns';
factor(idx.eta) = 1e6;
factor(idx.deltaC) = 1e9;

trueValue = state.true.*factor;
initialValue = state.initial.*factor;
estimatedValue = result.thetaHat.*factor;
errorValue = (result.thetaHat-state.true).*factor;
localStdBound = info.stdPhysicalSI.*factor;

tables.physical = table(labels, trueValue, initialValue, estimatedValue, ...
    errorValue, localStdBound, units, ...
    'VariableNames', {'Parameter', 'TrueValue', 'InitialValue', ...
    'EstimatedValue', 'EstimationError', 'LocalStdBound', 'Unit'});

path = (0:idx.K).';
tables.channel = table(path, real(cfg.true.beta(:)), imag(cfg.true.beta(:)), ...
    real(result.betaHat(:)), imag(result.betaHat(:)), ...
    abs(result.betaHat(:)-cfg.true.beta(:)), ...
    'VariableNames', {'Path', 'TrueReal', 'TrueImag', 'EstimatedReal', ...
    'EstimatedImag', 'AbsoluteError'});

uError = norm(result.thetaHat(idx.u0)-state.true(idx.u0));
vError = norm(result.thetaHat(idx.uDot)-state.true(idx.uDot));
etaErrorPpm = abs(result.thetaHat(idx.eta)-state.true(idx.eta))*1e6;
deltaErrorNs = abs(result.thetaHat(idx.deltaC)-state.true(idx.deltaC))*1e9;
normalizedDataCost = result.dataCost/numel(result.residual);

metric = {'UE position error'; 'UE velocity error'; ...
    'Absolute eta error'; 'Absolute delta_c error'; ...
    'Final objective'; 'Data objective'; 'Normalized data objective'; ...
    'Prior objective'; ...
    'Equivalent FIM rank'; 'Equivalent FIM dimension'; ...
    'Equivalent FIM condition number'};
value = [uError; vError; etaErrorPpm; deltaErrorNs; ...
    result.totalCost; result.dataCost; normalizedDataCost; result.priorCost; ...
    info.rank; info.dimension; info.conditionNumber];
unit = {'m'; 'm/s'; 'ppm'; 'ns'; '1'; '1'; '1/sample'; ...
    '1'; '1'; '1'; '1'};
tables.summary = table(metric, value, unit, ...
    'VariableNames', {'Metric', 'Value', 'Unit'});
end


function labels = physical_labels(idx)
labels = cell(idx.nPhysical, 1);
labels(idx.u0) = {'u0_x'; 'u0_y'; 'u0_z'};
labels(idx.uDot) = {'uDot_x'; 'uDot_y'; 'uDot_z'};
for k = 1:idx.K
    labels(idx.p0(:, k)) = {sprintf('p%d_x', k); ...
        sprintf('p%d_y', k); sprintf('p%d_z', k)};
    labels(idx.pDot(:, k)) = {sprintf('pDot%d_x', k); ...
        sprintf('pDot%d_y', k); sprintf('pDot%d_z', k)};
end
labels{idx.eta} = 'eta';
labels{idx.deltaC} = 'delta_c';
end
