function fig = plot_fim_diagnostics(results)
% Plot CRLB correlation and FIM diagnostics.

R = results.default;
S = results.sweep;
idx = results.indices;
K = size(results.config.risPositionM, 2);

labels = make_physical_labels(idx, K);
C = R.CphysicalSI;
stdVector = sqrt(max(diag(C), 0));
denominator = stdVector*stdVector.';
correlation = C ./ max(denominator, realmin);
correlation = max(min(correlation, 1), -1);
eigenvalues = sort(real(eig(R.JequivalentScaled)), 'descend');

fig = figure('Color', 'w', 'Name', 'FIM diagnostics', ...
    'Position', [70 70 1380 620]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(correlation);
axis image;
colormap(gca, parula);
colorbar;
caxis([-1, 1]);
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
    'XTickLabelRotation', 90, 'YTick', 1:numel(labels), ...
    'YTickLabel', labels, 'FontSize', 7);
title('CRLB correlation matrix');

nexttile;
semilogy(1:numel(eigenvalues), max(eigenvalues, realmin), 'o-', ...
    'LineWidth', 1.4);
grid on;
xlabel('Eigenvalue index'); ylabel('Scaled FIM eigenvalue');
title('Equivalent-FIM spectrum');

nexttile;
semilogy(S.snrDb, S.conditionEquivalent, 's-', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)'); ylabel('Condition number');
title('Numerical/identifiability condition');
title(tl, 'Information content and parameter coupling');
end


function labels = make_physical_labels(idx, K)
% Create physical-parameter labels.
labels = cell(idx.nPhysical, 1);
labels(idx.uePosition) = {'u_x'; 'u_y'; 'u_z'};
labels(idx.ueVelocity) = {'v_x'; 'v_y'; 'v_z'};

for k = 1:K
    labels(idx.risPosition(:, k)) = {sprintf('p%dx', k); ...
        sprintf('p%dy', k); sprintf('p%dz', k)};
    labels(idx.risVelocity(:, k)) = {sprintf('vp%dx', k); ...
        sprintf('vp%dy', k); sprintf('vp%dz', k)};
end

labels{idx.eta} = 'eta';
labels{idx.delta} = 'delta';
end
