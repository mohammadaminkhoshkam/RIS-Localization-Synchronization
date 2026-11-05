function figures = ris_plot_results(result, state, data, info, cfg, idx)
% Plot geometry, convergence, residuals, and channel estimates.

trueState = ris_unpack_state(state.true, idx);
initialState = ris_unpack_state(state.initial, idx);
estimatedState = ris_unpack_state(result.thetaHat, idx);

figures.geometry = figure('Color', 'w', 'Name', 'Estimator geometry');
hold on; grid on; axis equal; view(38, 24);
plot3(cfg.q(1), cfg.q(2), cfg.q(3), 'kp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'BS');
plot3(trueState.u0(1), trueState.u0(2), trueState.u0(3), 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'True UE');
plot3(initialState.u0(1), initialState.u0(2), initialState.u0(3), 'mo', ...
    'MarkerSize', 9, 'DisplayName', 'Initial UE');
plot3(estimatedState.u0(1), estimatedState.u0(2), estimatedState.u0(3), 'ro', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated UE');
plot3(trueState.p0(1, :), trueState.p0(2, :), trueState.p0(3, :), 'bs', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'b', 'DisplayName', 'True RIS');
plot3(cfg.prior.pMean(1, :), cfg.prior.pMean(2, :), cfg.prior.pMean(3, :), ...
    'cd', 'MarkerSize', 8, 'DisplayName', 'RIS prior mean');
plot3(estimatedState.p0(1, :), estimatedState.p0(2, :), ...
    estimatedState.p0(3, :), 'rs', 'MarkerSize', 9, ...
    'DisplayName', 'Estimated RIS');

plot3([cfg.q(1), trueState.u0(1)], [cfg.q(2), trueState.u0(2)], ...
    [cfg.q(3), trueState.u0(3)], 'k-', 'HandleVisibility', 'off');
for k = 1:idx.K
    plot3([cfg.q(1), trueState.p0(1, k), trueState.u0(1)], ...
        [cfg.q(2), trueState.p0(2, k), trueState.u0(2)], ...
        [cfg.q(3), trueState.p0(3, k), trueState.u0(3)], ...
        '--', 'HandleVisibility', 'off');
end
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('True, prior, initial, and estimated geometry');
legend('Location', 'bestoutside');

figures.diagnostics = figure('Color', 'w', 'Name', 'Estimator diagnostics', ...
    'Position', [100 100 1180 720]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
history = result.histories{result.bestStart};
if isempty(history)
    history = result.totalCost;
end
semilogy(0:numel(history)-1, max(history, realmin), 'o-', 'LineWidth', 1.4);
grid on; xlabel('Iteration'); ylabel('Objective');
title(sprintf('Best-start convergence (start %d)', result.bestStart));

nexttile;
residualMatrix = reshape(result.residual, cfg.N, cfg.L);
imagesc(data.times*1e3, data.n, abs(residualMatrix));
axis xy; colorbar; xlabel('Snapshot time (ms)'); ylabel('Subcarrier index');
title('|y-A(thetaHat)betaHat|');

nexttile;
plot(real(cfg.true.beta), imag(cfg.true.beta), 'bo', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'b'); hold on;
plot(real(result.betaHat), imag(result.betaHat), 'rx', ...
    'MarkerSize', 10, 'LineWidth', 1.8);
for k = 1:idx.M
    plot([real(cfg.true.beta(k)), real(result.betaHat(k))], ...
        [imag(cfg.true.beta(k)), imag(result.betaHat(k))], 'k:');
end
grid on; axis equal; xlabel('Real'); ylabel('Imaginary');
legend('True beta', 'Estimated beta', 'Location', 'best');
title('Complex path gains');

nexttile;
scaledError = abs(result.thetaHat-state.true)./state.scale;
bar(scaledError); grid on; xlabel('Physical parameter index');
ylabel('Absolute scaled error');
title(sprintf('FIM rank %d/%d, cond %.3g', ...
    info.rank, info.dimension, info.conditionNumber));
title(tl, 'Hybrid MAP/ML estimator diagnostics');

if result.syncInitialization.enabled
    s = result.syncInitialization;
    figures.initialization = figure('Color', 'w', ...
        'Name', 'Synchronization initialization', ...
        'Position', [130 130 1080 430]);
    tlSync = tiledlayout(1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    imagesc(s.coarseDeltaNs, s.coarseEtaPpm, ...
        log10(max(s.coarseCosts, realmin)));
    axis xy; colorbar; hold on;
    plot(s.bestDeltaNs, s.bestEtaPpm, 'rx', ...
        'MarkerSize', 11, 'LineWidth', 2);
    xlabel('delta_c (ns)'); ylabel('eta (ppm)');
    title('Coarse synchronization search');

    nexttile;
    imagesc(s.fineDeltaNs, s.fineEtaPpm, ...
        log10(max(s.fineCosts, realmin)));
    axis xy; colorbar; hold on;
    plot(s.bestDeltaNs, s.bestEtaPpm, 'rx', ...
        'MarkerSize', 11, 'LineWidth', 2);
    xlabel('delta_c (ns)'); ylabel('eta (ppm)');
    title('Fine synchronization search');
    title(tlSync, 'log10 concentrated MAP/ML objective');
end
end
