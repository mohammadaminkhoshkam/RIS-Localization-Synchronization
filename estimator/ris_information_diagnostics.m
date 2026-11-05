function info = ris_information_diagnostics(result, data, state, cfg, idx)
% Build the local hybrid FIM at the estimate.

Dall = [result.Dtheta, result.Ahat, 1j*result.Ahat];
scaleAll = [state.scale; ones(2*idx.M, 1)];
Dscaled = bsxfun(@times, Dall, scaleAll.');
J = (2/data.sigma2)*real(Dscaled'*Dscaled);

if cfg.prior.enabled
    for j = idx.p0(:).'
        J(j, j) = J(j, j) + (state.scale(j)/cfg.prior.positionStd)^2;
    end
    for j = idx.pDot(:).'
        J(j, j) = J(j, j) + (state.scale(j)/cfg.prior.velocityStd)^2;
    end
end

physical = 1:idx.nPhysical;
betaIndex = idx.nPhysical+(1:2*idx.M);
Jxx = J(physical, physical);
Jxb = J(physical, betaIndex);
Jbb = J(betaIndex, betaIndex);
Je = Jxx - Jxb*pinv(Jbb)*Jxb.';
Je = real((Je+Je.')/2);

[V, D] = eig(Je);
d = real(diag(D));
dMax = max(abs(d));
tol = max(size(Je))*eps(max(dMax, 1));
positive = d > tol;
rankJe = sum(positive);

invD = zeros(size(d));
invD(positive) = 1./d(positive);
Cscaled = V*diag(invD)*V.';
Csi = (state.scale*state.scale.').*Cscaled;
Csi = real((Csi+Csi.')/2);

info.JgammaScaled = J;
info.JequivalentScaled = Je;
info.CphysicalSI = Csi;
info.stdPhysicalSI = sqrt(max(diag(Csi), 0));
info.rank = rankJe;
info.dimension = idx.nPhysical;
info.eigenvalues = sort(d, 'descend');
if any(positive)
    info.conditionNumber = max(d(positive))/min(d(positive));
else
    info.conditionNumber = Inf;
end
info.fullRank = rankJe == idx.nPhysical;

if ~info.fullRank
    warning(['Equivalent FIM rank is %d of %d. Some estimates are not ' ...
        'locally identifiable; pseudo-inverse bounds are not finite CRLBs.'], ...
        info.rank, info.dimension);
end
end
