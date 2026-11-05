function data = ris_simulate_data(thetaTrue, betaTrue, cfg, idx)
% Generate pilots, noiseless samples, and complex Gaussian noise.

rng(cfg.randomSeed, 'twister');

if mod(cfg.N, 2) == 1
    data.n = (-(cfg.N-1)/2:(cfg.N-1)/2).';
else
    data.n = (-cfg.N/2:cfg.N/2-1).';
end

if cfg.centerTimes
    ell = (0:cfg.L-1) - (cfg.L-1)/2;
else
    ell = 0:cfg.L-1;
end
data.times = ell(:)*cfg.snapshotInterval;

qpskIndex = randi([0, 3], cfg.N, cfg.L);
data.pilots = exp(1j*(pi/4+(pi/2)*qpskIndex));

[data.Atrue, data.aTrue] = ris_observation_matrix( ...
    thetaTrue, data.pilots, data.n, data.times, cfg, idx);
data.mu = data.Atrue*betaTrue(:);
data.signalPower = mean(abs(data.mu).^2);
data.sigma2 = data.signalPower/10^(cfg.snrDb/10);
noise = sqrt(data.sigma2/2)*(randn(size(data.mu))+1j*randn(size(data.mu)));
data.y = data.mu + noise;
data.noise = noise;
end
