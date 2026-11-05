function fig = plot_system_geometry(cfg)
% Plot the BS, UE, RIS nodes, and paths.

K = size(cfg.risPositionM, 2);
q = cfg.bsPositionM(:);
u = cfg.uePositionM(:);
p = cfg.risPositionM;

fig = figure('Color', 'w', 'Name', 'System geometry');
hold on;
grid on;
axis equal;
view(38, 24);

plot3(q(1), q(2), q(3), 'kp', 'MarkerSize', 15, ...
    'MarkerFaceColor', [0.15 0.15 0.15], 'DisplayName', 'BS');
plot3(u(1), u(2), u(3), 'ro', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.85 0.20 0.20], 'DisplayName', 'UE');
plot3(p(1, :), p(2, :), p(3, :), 'bs', 'MarkerSize', 9, ...
    'MarkerFaceColor', [0.20 0.45 0.90], 'DisplayName', 'RIS');

plot3([q(1), u(1)], [q(2), u(2)], [q(3), u(3)], 'k-', ...
    'LineWidth', 1.8, 'DisplayName', 'LOS path');
for k = 1:K
    plot3([q(1), p(1, k), u(1)], [q(2), p(2, k), u(2)], ...
        [q(3), p(3, k), u(3)], '--', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('RIS path %d', k));
    text(p(1, k), p(2, k), p(3, k), sprintf('  RIS %d', k));
end

quiver3(u(1), u(2), u(3), cfg.ueVelocityMps(1), ...
    cfg.ueVelocityMps(2), cfg.ueVelocityMps(3), 2, ...
    'Color', [0.75 0 0], 'LineWidth', 1.7, 'MaxHeadSize', 0.8, ...
    'DisplayName', 'UE velocity');

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('BS, UE, RIS nodes, and propagation paths');
legend('Location', 'bestoutside');
end
