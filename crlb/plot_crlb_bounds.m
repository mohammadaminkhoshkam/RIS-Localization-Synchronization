function figures = plot_crlb_bounds(results)
% Plot aggregate and component CRLB results.

R = results.default;
S = results.sweep;
idx = results.indices;
K = size(results.config.risPositionM, 2);
M = K + 1;
xyz = {'x', 'y', 'z'};

%% Aggregate bounds
figures.mainBounds = figure('Color', 'w', ...
    'Name', 'Main bounds versus SNR', 'Position', [100 100 1150 760]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
semilogy(S.snrDb, S.uePEB_m, 'o-', 'LineWidth', 1.6); hold on;
semilogy(S.snrDb, S.ueVEB_mps, 's-', 'LineWidth', 1.6);
grid on;
xlabel('SNR (dB)'); ylabel('Lower bound');
legend('UE PEB (m)', 'UE VEB (m/s)', 'Location', 'southwest');
title('UE aggregate bounds');

nexttile;
for k = 1:K
    semilogy(S.snrDb, S.risPEB_m(k, :), 'o-', ...
        'LineWidth', 1.4); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('RIS PEB (m)');
legend(make_ris_labels(K), 'Location', 'best');
title('RIS position bounds');

nexttile;
for k = 1:K
    semilogy(S.snrDb, S.risVEB_mps(k, :), 's-', ...
        'LineWidth', 1.4); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('RIS VEB (m/s)');
legend(make_ris_labels(K), 'Location', 'best');
title('RIS velocity bounds');

nexttile;
semilogy(S.snrDb, S.etaStdPpm, 'd-', 'LineWidth', 1.6); hold on;
semilogy(S.snrDb, S.deltaStdNs, '^-', 'LineWidth', 1.6);
grid on;
xlabel('SNR (dB)'); ylabel('Std. lower bound');
legend('eta (ppm)', 'delta_c (ns)', 'Location', 'southwest');
title('Synchronization bounds');
title(tl, 'Aggregate hybrid CRLB results');

%% UE components
figures.ueComponents = figure('Color', 'w', ...
    'Name', 'UE component bounds', 'Position', [120 120 1050 690]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
for d = 1:3
    semilogy(S.snrDb, S.stdPhysicalSI(idx.uePosition(d), :), ...
        'LineWidth', 1.6); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('Position std. LB (m)');
legend(xyz, 'Location', 'southwest');
title('UE position components');

nexttile;
for d = 1:3
    semilogy(S.snrDb, S.stdPhysicalSI(idx.ueVelocity(d), :), ...
        'LineWidth', 1.6); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('Velocity std. LB (m/s)');
legend({'v_x', 'v_y', 'v_z'}, 'Location', 'southwest');
title('UE velocity components');
title(tl, 'Component-wise UE bounds');

%% RIS components
figures.risComponents = figure('Color', 'w', ...
    'Name', 'All RIS component bounds', 'Position', [80 40 1220 290*K]);
tl = tiledlayout(K, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:K
    nexttile;
    for d = 1:3
        semilogy(S.snrDb, S.stdPhysicalSI(idx.risPosition(d, k), :), ...
            'LineWidth', 1.35); hold on;
    end
    grid on;
    xlabel('SNR (dB)'); ylabel('Std. LB (m)');
    legend(xyz, 'Location', 'best');
    title(sprintf('RIS %d position components', k));

    nexttile;
    for d = 1:3
        semilogy(S.snrDb, S.stdPhysicalSI(idx.risVelocity(d, k), :), ...
            'LineWidth', 1.35); hold on;
    end
    grid on;
    xlabel('SNR (dB)'); ylabel('Std. LB (m/s)');
    legend({'v_x', 'v_y', 'v_z'}, 'Location', 'best');
    title(sprintf('RIS %d velocity components', k));
end
title(tl, 'All RIS scalar-parameter bounds');

%% Synchronization and channel
figures.synchronization = figure('Color', 'w', ...
    'Name', 'Synchronization and beta bounds', ...
    'Position', [90 90 1180 730]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
semilogy(S.snrDb, S.stdPhysicalSI(idx.eta, :)*1e6, 'd-', ...
    'LineWidth', 1.6);
grid on;
xlabel('SNR (dB)'); ylabel('Std. LB (ppm)');
title('Relative frequency error eta');

nexttile;
semilogy(S.snrDb, S.stdPhysicalSI(idx.delta, :)*1e9, '^-', ...
    'LineWidth', 1.6);
grid on;
xlabel('SNR (dB)'); ylabel('Std. LB (ns)');
title('Clock offset delta_c');

nexttile;
for k = 1:M
    semilogy(S.snrDb, S.stdAllSI(idx.betaReal(k), :), ...
        'LineWidth', 1.35); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('Std. LB');
legend(make_beta_labels(M, 'Re'), 'Location', 'best');
title('Real parts of beta_k');

nexttile;
for k = 1:M
    semilogy(S.snrDb, S.stdAllSI(idx.betaImag(k), :), ...
        'LineWidth', 1.35); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('Std. LB');
legend(make_beta_labels(M, 'Im'), 'Location', 'best');
title('Imaginary parts of beta_k');
title(tl, 'Synchronization and channel-coefficient bounds');

%% Default-SNR bars
figures.defaultBars = figure('Color', 'w', ...
    'Name', 'Default-SNR bounds', 'Position', [100 100 1180 720]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
bar(R.stdPhysicalSI(idx.uePosition));
grid on;
set(gca, 'XTick', 1:3, 'XTickLabel', xyz);
ylabel('Std. LB (m)');
title('UE position coordinates');

nexttile;
bar(R.stdPhysicalSI(idx.ueVelocity));
grid on;
set(gca, 'XTick', 1:3, 'XTickLabel', {'v_x', 'v_y', 'v_z'});
ylabel('Std. LB (m/s)');
title('UE velocity coordinates');

nexttile;
bar(1:K, [R.risPEB_m(:), R.risVEB_mps(:)], 'grouped');
grid on;
xlabel('RIS index'); ylabel('Aggregate lower bound');
legend('PEB (m)', 'VEB (m/s)', 'Location', 'best');
title('RIS aggregate bounds');

nexttile;
bar(0:M-1, [R.betaRealStd(:), R.betaImagStd(:)], 'grouped');
grid on;
xlabel('Path index k'); ylabel('Std. LB');
legend('Re{beta_k}', 'Im{beta_k}', 'Location', 'best');
title('Channel coefficients');
title(tl, sprintf('Bounds at SNR = %.1f dB', results.config.snrDb));
end


function labels = make_ris_labels(K)
% Create RIS legend labels.
labels = cell(1, K);
for k = 1:K
    labels{k} = sprintf('RIS %d', k);
end
end


function labels = make_beta_labels(M, partName)
% Create channel-coefficient labels.
labels = cell(1, M);
for k = 1:M
    labels{k} = sprintf('%s{beta_%d}', partName, k-1);
end
end
