function s = ris_unpack_state(theta, idx)
% Unpack the physical parameter vector.

s.u0 = theta(idx.u0);
s.uDot = theta(idx.uDot);
s.p0 = reshape(theta(idx.p0(:)), 3, idx.K);
s.pDot = reshape(theta(idx.pDot(:)), 3, idx.K);
s.eta = theta(idx.eta);
s.deltaC = theta(idx.deltaC);
end
