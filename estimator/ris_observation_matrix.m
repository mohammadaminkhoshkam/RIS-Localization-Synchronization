function [A, a] = ris_observation_matrix(theta, pilots, n, times, cfg, idx)
% Build the OFDM observation matrix A(theta).

a = ris_effective_delays(theta, times, cfg, idx);
A = zeros(cfg.N*cfg.L, idx.M);
gain = sqrt(cfg.energy/cfg.N);

for k = 1:idx.M
    phase = exp(-1j*2*pi*cfg.deltaF*(n*a(:, k).'));
    Ak = gain*pilots.*phase;
    A(:, k) = Ak(:);
end
end
