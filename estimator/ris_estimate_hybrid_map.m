function result = ris_estimate_hybrid_map(data, state, cfg, idx)
% Run bounded multi-start variable projection.

scale = state.scale;
z0 = state.initial./scale;
[searchLower, searchUpper] = make_search_bounds(state, cfg, idx);
zLower = searchLower./scale;
zUpper = searchUpper./scale;
startStdZ = state.startStd./scale;

objective = @(z) ris_hybrid_map_cost(z, data.y, data.pilots, ...
    data.n, data.times, data.sigma2, cfg, idx, scale);

zBase = clamp(z0, zLower, zUpper);
[zBase, syncInfo] = ris_coarse_sync_search( ...
    zBase, objective, zLower, zUpper, cfg, idx, scale);

rng(cfg.randomSeed + 1, 'twister');
nStarts = cfg.optimizer.numStarts;
nCandidates = max(nStarts, ...
    ceil(cfg.optimizer.candidateMultiplier*nStarts));
candidatePoints = zeros(idx.nPhysical, nCandidates);
candidateCosts = Inf(nCandidates, 1);
candidatePoints(:, 1) = zBase;
candidateCosts(1) = objective(zBase);
for c = 2:nCandidates
    trial = zBase + startStdZ.*randn(size(zBase));
    candidatePoints(:, c) = clamp(trial, zLower, zUpper);
    candidateCosts(c) = objective(candidatePoints(:, c));
end

[~, candidateOrder] = sort(candidateCosts, 'ascend');
candidateOrder(candidateOrder == 1) = [];
selected = 1;
if nStarts > 1
    selected = [1; candidateOrder(1:nStarts-1)];
end
startPoints = candidatePoints(:, selected);
startCosts = candidateCosts(selected);

bestCost = Inf;
bestZ = startPoints(:, 1);
bestExitFlag = NaN;
bestOutput = struct();
bestStart = 1;
histories = cell(nStarts, 1);
finalCosts = Inf(nStarts, 1);
exitFlags = NaN(nStarts, 1);

hasFmincon = exist('fmincon', 'file') == 2;
for s = 1:nStarts
    history = [];
    if hasFmincon
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...
            'SpecifyObjectiveGradient', true, ...
            'MaxIterations', cfg.optimizer.maxIterations, ...
            'MaxFunctionEvaluations', cfg.optimizer.maxFunctionEvaluations, ...
            'OptimalityTolerance', cfg.optimizer.optimalityTolerance, ...
            'StepTolerance', cfg.optimizer.stepTolerance, ...
            'FunctionTolerance', 1e-10, ...
            'ConstraintTolerance', 1e-10, ...
            'TypicalX', max(abs(zBase), 1), ...
            'Display', cfg.optimizer.display, ...
            'OutputFcn', @record_history);
        [zHat, fval, exitFlag, output] = fmincon( ...
            objective, startPoints(:, s), [], [], [], [], ...
            zLower, zUpper, [], options);
    else
        warning('fmincon was not found. Using penalized fminsearch.');
        penalty = @(z) bounded_objective(z, objective, zLower, zUpper);
        options = optimset('MaxIter', cfg.optimizer.maxIterations, ...
            'MaxFunEvals', cfg.optimizer.maxFunctionEvaluations, ...
            'TolX', cfg.optimizer.stepTolerance, ...
            'TolFun', cfg.optimizer.optimalityTolerance, ...
            'Display', cfg.optimizer.display);
        [zHat, fval, exitFlag, output] = fminsearch( ...
            penalty, startPoints(:, s), options);
        zHat = clamp(zHat, zLower, zUpper);
        fval = objective(zHat);
        history = fval;
    end

    finalCosts(s) = fval;
    exitFlags(s) = exitFlag;
    histories{s} = history;
    if fval < bestCost
        bestCost = fval;
        bestZ = zHat;
        bestExitFlag = exitFlag;
        bestOutput = output;
        bestStart = s;
    end
end

[~, ~, details] = objective(bestZ);
result.thetaHat = details.theta;
result.betaHat = details.beta;
result.residual = details.residual;
result.Ahat = details.A;
result.aHat = details.a;
result.Dtheta = details.Dtheta;
result.totalCost = bestCost;
result.dataCost = details.dataCost;
result.priorCost = details.priorCost;
result.bestStart = bestStart;
result.exitFlag = bestExitFlag;
result.output = bestOutput;
result.histories = histories;
result.startPoints = bsxfun(@times, startPoints, scale);
result.startCosts = startCosts;
result.finalCosts = finalCosts;
result.exitFlags = exitFlags;
result.usedFmincon = hasFmincon;
result.syncInitialization = syncInfo;
result.searchLower = searchLower;
result.searchUpper = searchUpper;

    function stop = record_history(~, optimValues, stateName)
        stop = false;
        if strcmp(stateName, 'init') || strcmp(stateName, 'iter')
            history(end+1, 1) = optimValues.fval;
        end
    end
end


function [lower, upper] = make_search_bounds(state, cfg, idx)
% Restrict the local physical search region.

lower = state.lower;
upper = state.upper;
if ~cfg.optimizer.localBounds.enabled
    return;
end

u0Radius = expand_radius(cfg.optimizer.localBounds.u0Radius, 3);
uDotRadius = expand_radius(cfg.optimizer.localBounds.uDotRadius, 3);
p0Radius = expand_radius(cfg.optimizer.localBounds.p0Radius, numel(idx.p0));
pDotRadius = expand_radius( ...
    cfg.optimizer.localBounds.pDotRadius, numel(idx.pDot));

lower(idx.u0) = max(lower(idx.u0), state.initial(idx.u0)-u0Radius);
upper(idx.u0) = min(upper(idx.u0), state.initial(idx.u0)+u0Radius);
lower(idx.uDot) = max( ...
    lower(idx.uDot), state.initial(idx.uDot)-uDotRadius);
upper(idx.uDot) = min( ...
    upper(idx.uDot), state.initial(idx.uDot)+uDotRadius);
lower(idx.p0(:)) = max( ...
    lower(idx.p0(:)), state.initial(idx.p0(:))-p0Radius);
upper(idx.p0(:)) = min( ...
    upper(idx.p0(:)), state.initial(idx.p0(:))+p0Radius);
lower(idx.pDot(:)) = max( ...
    lower(idx.pDot(:)), state.initial(idx.pDot(:))-pDotRadius);
upper(idx.pDot(:)) = min( ...
    upper(idx.pDot(:)), state.initial(idx.pDot(:))+pDotRadius);
end


function radius = expand_radius(value, count)
% Expand a scalar or validate a radius vector.

if isscalar(value)
    radius = value*ones(count, 1);
else
    radius = value(:);
    if numel(radius) ~= count
        error('A local-bound radius has an invalid size.');
    end
end
end


function value = bounded_objective(z, objective, lower, upper)
zClamped = clamp(z, lower, upper);
violation = z-zClamped;
value = objective(zClamped) + 1e8*sum(violation.^2);
end


function x = clamp(x, lower, upper)
x = min(max(x, lower), upper);
end
