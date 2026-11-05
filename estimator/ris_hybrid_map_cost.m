function [cost, gradZ, details] = ris_hybrid_map_cost( ...
        z, y, pilots, n, times, sigma2, cfg, idx, scale)
% Evaluate the concentrated hybrid MAP/ML objective.

theta = scale.*z;
[A, a] = ris_observation_matrix(theta, pilots, n, times, cfg, idx);
beta = stable_least_squares(A, y);
residual = y - A*beta;

dataCost = real(residual'*residual)/sigma2;
priorCost = 0;
priorGradient = zeros(idx.nPhysical, 1);

if cfg.prior.enabled
    pError = theta(idx.p0(:)) - cfg.prior.pMean(:);
    vError = theta(idx.pDot(:)) - cfg.prior.pDotMean(:);
    priorCost = 0.5*sum((pError/cfg.prior.positionStd).^2) ...
        + 0.5*sum((vError/cfg.prior.velocityStd).^2);
    priorGradient(idx.p0(:)) = pError/(cfg.prior.positionStd^2);
    priorGradient(idx.pDot(:)) = vError/(cfg.prior.velocityStd^2);
end

cost = dataCost + priorCost;

if nargout >= 2
    Dtheta = ris_mean_jacobian(theta, beta, A, a, ...
        pilots, n, times, cfg, idx);
    dataGradient = -(2/sigma2)*real(Dtheta'*residual);
    gradTheta = dataGradient + priorGradient;
    gradZ = scale.*gradTheta;
else
    Dtheta = [];
    gradZ = [];
end

if nargout >= 3
    details.theta = theta;
    details.beta = beta;
    details.A = A;
    details.a = a;
    details.residual = residual;
    details.dataCost = dataCost;
    details.priorCost = priorCost;
    details.Dtheta = Dtheta;
end
end


function beta = stable_least_squares(A, y)
[U, S, V] = svd(A, 'econ');
s = diag(S);
tol = max(size(A))*eps(max(s));
keep = s > tol;
if any(keep)
    beta = V(:, keep)*((U(:, keep)'*y)./s(keep));
else
    beta = zeros(size(A, 2), 1);
end
end
