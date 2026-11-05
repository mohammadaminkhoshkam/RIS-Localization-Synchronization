function [zBest, info] = ris_coarse_sync_search( ...
        z0, objective, zLower, zUpper, cfg, idx, scale)
% Refine eta and delta_c before the full optimization.

info.enabled = cfg.optimizer.syncSearch.enabled;
info.initialCost = objective(z0);
zBest = z0;

if ~info.enabled
    info.bestCost = info.initialCost;
    info.initialEtaPpm = scale(idx.eta)*z0(idx.eta)*1e6;
    info.initialDeltaNs = scale(idx.deltaC)*z0(idx.deltaC)*1e9;
    info.bestEtaPpm = info.initialEtaPpm;
    info.bestDeltaNs = info.initialDeltaNs;
    return;
end

s = cfg.optimizer.syncSearch;
etaHalf = s.etaHalfSpan/scale(idx.eta);
deltaHalf = s.deltaHalfSpan/scale(idx.deltaC);

etaCoarse = linspace(max(zLower(idx.eta), z0(idx.eta)-etaHalf), ...
    min(zUpper(idx.eta), z0(idx.eta)+etaHalf), s.coarseEtaPoints);
deltaCoarse = linspace(max(zLower(idx.deltaC), ...
    z0(idx.deltaC)-deltaHalf), min(zUpper(idx.deltaC), ...
    z0(idx.deltaC)+deltaHalf), s.coarseDeltaPoints);

[zBest, coarseCosts] = search_grid( ...
    z0, objective, etaCoarse, deltaCoarse, idx);

etaStep = grid_step(etaCoarse, etaHalf);
deltaStep = grid_step(deltaCoarse, deltaHalf);
etaFine = linspace(max(zLower(idx.eta), zBest(idx.eta)-etaStep), ...
    min(zUpper(idx.eta), zBest(idx.eta)+etaStep), s.fineEtaPoints);
deltaFine = linspace(max(zLower(idx.deltaC), ...
    zBest(idx.deltaC)-deltaStep), min(zUpper(idx.deltaC), ...
    zBest(idx.deltaC)+deltaStep), s.fineDeltaPoints);

[zBest, fineCosts] = search_grid( ...
    zBest, objective, etaFine, deltaFine, idx);

info.bestCost = objective(zBest);
info.initialEtaPpm = scale(idx.eta)*z0(idx.eta)*1e6;
info.initialDeltaNs = scale(idx.deltaC)*z0(idx.deltaC)*1e9;
info.bestEtaPpm = scale(idx.eta)*zBest(idx.eta)*1e6;
info.bestDeltaNs = scale(idx.deltaC)*zBest(idx.deltaC)*1e9;
info.coarseEtaPpm = scale(idx.eta)*etaCoarse*1e6;
info.coarseDeltaNs = scale(idx.deltaC)*deltaCoarse*1e9;
info.coarseCosts = coarseCosts;
info.fineEtaPpm = scale(idx.eta)*etaFine*1e6;
info.fineDeltaNs = scale(idx.deltaC)*deltaFine*1e9;
info.fineCosts = fineCosts;
end


function [zBest, costs] = search_grid( ...
        zCenter, objective, etaGrid, deltaGrid, idx)
% Evaluate one eta-delta grid.

costs = Inf(numel(etaGrid), numel(deltaGrid));
bestCost = Inf;
zBest = zCenter;
for i = 1:numel(etaGrid)
    for j = 1:numel(deltaGrid)
        zTrial = zCenter;
        zTrial(idx.eta) = etaGrid(i);
        zTrial(idx.deltaC) = deltaGrid(j);
        costs(i, j) = objective(zTrial);
        if costs(i, j) < bestCost
            bestCost = costs(i, j);
            zBest = zTrial;
        end
    end
end
end


function step = grid_step(grid, fallback)
% Return a positive grid spacing.

if numel(grid) > 1
    step = max(abs(diff(grid)));
else
    step = fallback;
end
end
