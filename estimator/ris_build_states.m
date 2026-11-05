function state = ris_build_states(cfg, idx)
% Pack true, initial, and bounded states.

state.true = pack_state(cfg.true, idx);
state.initial = pack_state(cfg.initial, idx);

state.lower = zeros(idx.nPhysical, 1);
state.upper = zeros(idx.nPhysical, 1);
state.lower(idx.u0) = cfg.bounds.u0Lower;
state.upper(idx.u0) = cfg.bounds.u0Upper;
state.lower(idx.uDot) = cfg.bounds.uDotLower;
state.upper(idx.uDot) = cfg.bounds.uDotUpper;
state.lower(idx.p0(:)) = cfg.bounds.p0Lower(:);
state.upper(idx.p0(:)) = cfg.bounds.p0Upper(:);
state.lower(idx.pDot(:)) = cfg.bounds.pDotLower(:);
state.upper(idx.pDot(:)) = cfg.bounds.pDotUpper(:);
state.lower(idx.eta) = cfg.bounds.etaLower;
state.upper(idx.eta) = cfg.bounds.etaUpper;
state.lower(idx.deltaC) = cfg.bounds.deltaLower;
state.upper(idx.deltaC) = cfg.bounds.deltaUpper;

state.scale = ones(idx.nPhysical, 1);
state.scale(idx.eta) = 1e-6;
state.scale(idx.deltaC) = 1e-9;

state.startStd = zeros(idx.nPhysical, 1);
state.startStd(idx.u0) = cfg.optimizer.startStd.u0;
state.startStd(idx.uDot) = cfg.optimizer.startStd.uDot;
state.startStd(idx.p0(:)) = cfg.optimizer.startStd.p0;
state.startStd(idx.pDot(:)) = cfg.optimizer.startStd.pDot;
state.startStd(idx.eta) = cfg.optimizer.startStd.eta;
state.startStd(idx.deltaC) = cfg.optimizer.startStd.deltaC;

if numel(cfg.true.beta) ~= idx.M
    error('cfg.true.beta must contain K+1=%d coefficients.', idx.M);
end
end


function theta = pack_state(s, idx)
theta = zeros(idx.nPhysical, 1);
theta(idx.u0) = s.u0(:);
theta(idx.uDot) = s.uDot(:);
theta(idx.p0(:)) = s.p0(:);
theta(idx.pDot(:)) = s.pDot(:);
theta(idx.eta) = s.eta;
theta(idx.deltaC) = s.deltaC;
end
