function [sensitivityTable, fig] = analyze_ris_prior(cfg)
% Evaluate CRLB sensitivity to RIS prior accuracy.

positionStdM = [0.02, 0.05, 0.10, 0.20, 0.50];
velocityStdMps = positionStdM/2;
uePEB = zeros(size(positionStdM));
ueVEB = zeros(size(positionStdM));
meanRisPEB = zeros(size(positionStdM));
meanRisVEB = zeros(size(positionStdM));

for i = 1:numel(positionStdM)
    cfgPrior = cfg;
    cfgPrior.risPrior.positionStdM = positionStdM(i);
    cfgPrior.risPrior.velocityStdMps = velocityStdMps(i);
    cfgPrior.snrSweepDb = cfgPrior.snrDb;
    cfgPrior.makePlots = false;
    cfgPrior.verbose = false;
    cfgPrior.saveResults = false;

    priorResults = simulate_all_crlb_bounds(cfgPrior);
    uePEB(i) = priorResults.default.uePEB_m;
    ueVEB(i) = priorResults.default.ueVEB_mps;
    meanRisPEB(i) = mean(priorResults.default.risPEB_m);
    meanRisVEB(i) = mean(priorResults.default.risVEB_mps);
end

sensitivityTable = table(positionStdM(:), velocityStdMps(:), ...
    uePEB(:), ueVEB(:), meanRisPEB(:), meanRisVEB(:), ...
    'VariableNames', {'RIS_PositionPriorStd_m', ...
    'RIS_VelocityPriorStd_mps', 'UE_PEB_m', 'UE_VEB_mps', ...
    'Mean_RIS_PEB_m', 'Mean_RIS_VEB_mps'});

fig = figure('Color', 'w', 'Name', 'RIS prior sensitivity', ...
    'Position', [130 130 1050 430]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
loglog(positionStdM, uePEB, 'o-', 'LineWidth', 1.6); hold on;
loglog(positionStdM, meanRisPEB, 's-', 'LineWidth', 1.6);
grid on;
xlabel('RIS position prior std. (m)'); ylabel('Position bound (m)');
legend('UE PEB', 'Mean RIS PEB', 'Location', 'northwest');
title('Position sensitivity to RIS prior');

nexttile;
loglog(velocityStdMps, ueVEB, 'o-', 'LineWidth', 1.6); hold on;
loglog(velocityStdMps, meanRisVEB, 's-', 'LineWidth', 1.6);
grid on;
xlabel('RIS velocity prior std. (m/s)');
ylabel('Velocity bound (m/s)');
legend('UE VEB', 'Mean RIS VEB', 'Location', 'northwest');
title('Velocity sensitivity to RIS prior');
title(tl, 'Effect of prior information in the hybrid FIM');
end
