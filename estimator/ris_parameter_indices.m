function idx = ris_parameter_indices(K)
% Index the physical parameter vector.

j = 1;
idx.u0 = j:j+2; j = j + 3;
idx.uDot = j:j+2; j = j + 3;
idx.p0 = reshape(j:j+3*K-1, 3, K); j = j + 3*K;
idx.pDot = reshape(j:j+3*K-1, 3, K); j = j + 3*K;
idx.eta = j; j = j + 1;
idx.deltaC = j; j = j + 1;
idx.nPhysical = j - 1;
idx.K = K;
idx.M = K + 1;
end
