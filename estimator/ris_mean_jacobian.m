function D = ris_mean_jacobian(theta, beta, A, a, pilots, n, times, cfg, idx)
% Differentiate A(theta)beta with respect to physical parameters.

D = zeros(cfg.N*cfg.L, idx.nPhysical);
h = zeros(idx.nPhysical, 1);
h(idx.u0) = cfg.finiteDifference.positionStep;
h(idx.uDot) = cfg.finiteDifference.velocityStep;
h(idx.p0(:)) = cfg.finiteDifference.positionStep;
h(idx.pDot(:)) = cfg.finiteDifference.velocityStep;

for j = 1:idx.nPhysical
    if j == idx.eta
        da = repmat(times(:), 1, idx.M);
    elseif j == idx.deltaC
        da = -ones(cfg.L, idx.M);
    else
        thetaPlus = theta;
        thetaMinus = theta;
        thetaPlus(j) = thetaPlus(j) + h(j);
        thetaMinus(j) = thetaMinus(j) - h(j);
        da = (ris_effective_delays(thetaPlus, times, cfg, idx) ...
            - ris_effective_delays(thetaMinus, times, cfg, idx))/(2*h(j));
    end

    dmu = zeros(cfg.N, cfg.L);
    for k = 1:idx.M
        Ak = reshape(A(:, k), cfg.N, cfg.L);
        phaseSlope = -1j*2*pi*cfg.deltaF*(n*da(:, k).');
        dmu = dmu + beta(k)*Ak.*phaseSlope;
    end
    D(:, j) = dmu(:);
end
end
