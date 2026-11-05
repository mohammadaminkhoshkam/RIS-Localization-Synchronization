function results = simulate_all_crlb_bounds(userCfg)

if nargin < 1 || isempty(userCfg)
    cfg = default_config();
else
    cfg = merge_struct(default_config(), userCfg);
end

validate_config(cfg);
rng(cfg.randomSeed, 'twister');

K = size(cfg.risPositionM, 2);       % Number of RIS paths
M = K + 1;                          % LOS plus K RIS paths
idx = make_indices(K);

% Balanced subcarrier indices. An odd N gives exact symmetry about zero.
if mod(cfg.numSubcarriers, 2) == 1
    n = (-(cfg.numSubcarriers - 1)/2:(cfg.numSubcarriers - 1)/2).';
else
    n = (-cfg.numSubcarriers/2:cfg.numSubcarriers/2 - 1).';
end

if cfg.centerObservationTimes
    ell = (0:cfg.numSnapshots - 1) - (cfg.numSnapshots - 1)/2;
else
    ell = 0:cfg.numSnapshots - 1;
end
timesSec = ell(:) * cfg.snapshotIntervalSec;

% Unit-modulus known QPSK pilots; Communications Toolbox is not required.
qpskIndex = randi([0, 3], cfg.numSubcarriers, cfg.numSnapshots);
pilots = exp(1j * (pi/4 + (pi/2) * qpskIndex));

phys = zeros(idx.nPhysical, 1);
phys(idx.uePosition) = cfg.uePositionM(:);
phys(idx.ueVelocity) = cfg.ueVelocityMps(:);
phys(idx.risPosition(:)) = cfg.risPositionM(:);
phys(idx.risVelocity(:)) = cfg.risVelocityMps(:);
phys(idx.eta) = cfg.eta;
phys(idx.delta) = cfg.clockOffsetSec;

beta = cfg.beta(:);
if numel(beta) ~= M
    error('cfg.beta must contain K+1=%d complex path coefficients.', M);
end

[mu, Dphysical, phase, aSec] = mean_and_jacobian( ...
    cfg, phys, beta, pilots, n, timesSec, idx);

% Append the analytic derivatives with respect to Re{beta_k}, Im{beta_k}.
D = zeros(numel(mu), idx.nAll);
D(:, 1:idx.nPhysical) = Dphysical;
for k = 1:M
    base = sqrt(cfg.ofdmEnergy/cfg.numSubcarriers) .* pilots .* phase(:, :, k);
    D(:, idx.betaReal(k)) = base(:);
    D(:, idx.betaImag(k)) = 1j * base(:);
end

% Work in well-scaled coordinates to avoid numerical problems caused by
% mixing metres, ppm, nanoseconds, and dimensionless channel coefficients.
% If gamma_SI = diag(scale) gamma_scaled, then D_scaled=D_SI*diag(scale).
scale = ones(idx.nAll, 1);
scale(idx.eta) = 1e-6;       % one scaled unit = 1 ppm
scale(idx.delta) = 1e-9;     % one scaled unit = 1 ns
Dscaled = bsxfun(@times, D, scale.');

% Data FIM before division by noise variance.
JdataUnitNoise = 2 * real(Dscaled' * Dscaled);
JdataUnitNoise = symmetrize(JdataUnitNoise);
JpriorScaled = make_prior_fim(cfg, scale, idx);

signalPower = mean(abs(mu(:)).^2);
if signalPower <= 0 || ~isfinite(signalPower)
    error('The default received signal has invalid or zero power.');
end

% Default-SNR result.
defaultBound = evaluate_at_snr(cfg.snrDb, signalPower, ...
    JdataUnitNoise, JpriorScaled, scale, idx, K);

% SNR sweep. The signal Jacobian is unchanged; only sigma^2 changes.
snrList = cfg.snrSweepDb(:).';
nSNR = numel(snrList);
sweep = initialize_sweep(K, M, nSNR, idx.nPhysical, idx.nAll);
for s = 1:nSNR
    b = evaluate_at_snr(snrList(s), signalPower, ...
        JdataUnitNoise, JpriorScaled, scale, idx, K);
    sweep = store_sweep_point(sweep, b, s);
end
sweep.snrDb = snrList;

parameterTable = make_parameter_table(defaultBound, idx, K);
sweepTable = make_sweep_table(sweep, K, M);

if cfg.makePlots
    [fig1, fig2, fig3] = make_plots(cfg, defaultBound, sweep, idx, K, M);
    figureHandles = [fig1, fig2, fig3];
else
    fig1 = [];
    fig2 = [];
    fig3 = [];
    figureHandles = [];
end

results = struct();
results.config = cfg;
results.indices = idx;
results.subcarrierIndices = n;
results.observationTimesSec = timesSec;
results.effectiveDelaySec = aSec;
results.noiselessMean = mu;
results.default = defaultBound;
results.sweep = sweep;
results.parameterTable = parameterTable;
results.sweepTable = sweepTable;
results.figureHandles = figureHandles;

if cfg.verbose
    print_summary(defaultBound, cfg, K);
end

if cfg.saveResults
    save_outputs(results, parameterTable, sweepTable, fig1, fig2, fig3, cfg);
end
end


function cfg = default_config()
% Numerically useful, reproducible defaults. All quantities use SI units.
cfg.speedOfLightMps = 299792458;
cfg.numSubcarriers = 129;             % odd -> symmetric n=-64:64
cfg.numSnapshots = 32;
cfg.subcarrierSpacingHz = 480e3;
cfg.snapshotIntervalSec = 5e-3;
cfg.centerObservationTimes = true;
cfg.ofdmEnergy = 1;
cfg.snrDb = 20;
cfg.snrSweepDb = -10:5:30;
cfg.randomSeed = 42;

% Geometry: q is the static BS. Columns of risPositionM are p_k.
cfg.bsPositionM = [0; 0; 5];
cfg.uePositionM = [8; 7; 1.5];
cfg.ueVelocityMps = [1.2; -0.4; 0.1];
cfg.risPositionM = [0, 20, 18; ...
                    18, 0, 18; ...
                    3, 4, 6];
cfg.risVelocityMps = zeros(3, 3);

cfg.eta = 2e-6;                       % 2 ppm relative frequency error
cfg.clockOffsetSec = 40e-9;           % 40 ns
cfg.beta = [1.00*exp(1j*0.15); ...    % beta_0: LOS
            0.65*exp(-1j*0.80); ...
            0.50*exp(1j*1.40); ...
            0.40*exp(-1j*2.20)];

% Prior precision in the PDF is Q. Here Q is formed from intuitive standard
% deviations. Set enabled=false to remove the RIS prior entirely.
cfg.risPrior.enabled = true;
cfg.risPrior.positionStdM = 0.10;
cfg.risPrior.velocityStdMps = 0.05;

cfg.geometryModel = 'reference-linearized';

% Central finite differences are applied to a_{l,k}; eta and delta_c use
% exact analytic derivatives. These steps are deliberately unit-aware.
cfg.finiteDifference.positionStepM = 1e-3;
cfg.finiteDifference.velocityStepMps = 1e-2;

cfg.makePlots = true;
cfg.verbose = true;
cfg.saveResults = true;
cfg.outputFolder = 'crlb_results';
end


function validate_config(cfg)
if cfg.numSubcarriers < 3 || cfg.numSnapshots < 2
    error('At least 3 subcarriers and 2 snapshots are required.');
end
if cfg.subcarrierSpacingHz <= 0 || cfg.snapshotIntervalSec <= 0
    error('Subcarrier spacing and snapshot interval must be positive.');
end
if size(cfg.risPositionM, 1) ~= 3
    error('cfg.risPositionM must be a 3-by-K matrix.');
end
if ~isequal(size(cfg.risVelocityMps), size(cfg.risPositionM))
    error('cfg.risVelocityMps must have the same size as cfg.risPositionM.');
end
if numel(cfg.bsPositionM) ~= 3 || numel(cfg.uePositionM) ~= 3 || ...
        numel(cfg.ueVelocityMps) ~= 3
    error('BS/UE position and UE velocity must be three-dimensional vectors.');
end
validModels = {'reference-linearized', 'time-varying'};
if ~any(strcmpi(cfg.geometryModel, validModels))
    error('cfg.geometryModel must be reference-linearized or time-varying.');
end
if cfg.risPrior.enabled && (cfg.risPrior.positionStdM <= 0 || ...
        cfg.risPrior.velocityStdMps <= 0)
    error('Enabled prior standard deviations must be positive.');
end
end


function idx = make_indices(K)
M = K + 1;
j = 1;
idx.uePosition = j:j+2; j = j + 3;
idx.ueVelocity = j:j+2; j = j + 3;
idx.risPosition = reshape(j:j+3*K-1, 3, K); j = j + 3*K;
idx.risVelocity = reshape(j:j+3*K-1, 3, K); j = j + 3*K;
idx.eta = j; j = j + 1;
idx.delta = j; j = j + 1;
idx.nPhysical = j - 1;
idx.betaReal = j:j+M-1; j = j + M;
idx.betaImag = j:j+M-1; j = j + M;
idx.nAll = j - 1;
end


function [mu, Dphysical, phase, aSec] = mean_and_jacobian( ...
        cfg, phys, beta, pilots, n, timesSec, idx)
N = cfg.numSubcarriers;
L = cfg.numSnapshots;
M = numel(beta);

aSec = effective_delays(cfg, phys, timesSec, idx);
phase = zeros(N, L, M);
for k = 1:M
    phase(:, :, k) = exp(-1j * 2*pi * cfg.subcarrierSpacingHz ...
        * (n * aSec(:, k).'));
end

pathSum = zeros(N, L);
for k = 1:M
    pathSum = pathSum + beta(k) * phase(:, :, k);
end
muMatrix = sqrt(cfg.ofdmEnergy/N) .* pilots .* pathSum;
mu = muMatrix(:);

Dphysical = zeros(N*L, idx.nPhysical);
h = zeros(idx.nPhysical, 1);
h(idx.uePosition) = cfg.finiteDifference.positionStepM;
h(idx.ueVelocity) = cfg.finiteDifference.velocityStepMps;
h(idx.risPosition(:)) = cfg.finiteDifference.positionStepM;
h(idx.risVelocity(:)) = cfg.finiteDifference.velocityStepMps;

for j = 1:idx.nPhysical
    if j == idx.eta
        % From a_lk=tau_lk-(nu_lk-eta)T_l-delta_c.
        da = repmat(timesSec(:), 1, M);
    elseif j == idx.delta
        da = -ones(L, M);
    else
        hp = h(j);
        physPlus = phys;
        physMinus = phys;
        physPlus(j) = physPlus(j) + hp;
        physMinus(j) = physMinus(j) - hp;
        da = (effective_delays(cfg, physPlus, timesSec, idx) ...
            - effective_delays(cfg, physMinus, timesSec, idx)) / (2*hp);
    end

    dmu = zeros(N, L);
    for k = 1:M
        phaseDerivative = phase(:, :, k) .* ...
            (-1j * 2*pi * cfg.subcarrierSpacingHz * (n * da(:, k).'));
        dmu = dmu + beta(k) * phaseDerivative;
    end
    dmu = sqrt(cfg.ofdmEnergy/N) .* pilots .* dmu;
    Dphysical(:, j) = dmu(:);
end
end


function aSec = effective_delays(cfg, phys, timesSec, idx)
% Equations for LOS (k=0) and RIS paths (k>=1) from page 5.
c = cfg.speedOfLightMps;
q = cfg.bsPositionM(:);
u0 = phys(idx.uePosition);
uDot = phys(idx.ueVelocity);
p0 = reshape(phys(idx.risPosition(:)), 3, []);
pDot = reshape(phys(idx.risVelocity(:)), 3, []);
eta = phys(idx.eta);
delta = phys(idx.delta);
K = size(p0, 2);
M = K + 1;
L = numel(timesSec);
aSec = zeros(L, M);

for ell = 1:L
    T = timesSec(ell);
    if strcmpi(cfg.geometryModel, 'time-varying')
        u = u0 + uDot*T;
        p = p0 + pDot*T;
    else
        u = u0;
        p = p0;
    end

    % k=0, LOS path.
    d0 = u - q;
    r0 = norm(d0);
    check_distance(r0, 'UE and BS');
    tau0 = r0/c;
    nu0 = (uDot.' * d0)/(c*r0);
    aSec(ell, 1) = tau0 - (nu0 - eta)*T - delta;

    % k>=1, BS -> RIS_k -> UE.
    for k = 1:K
        dBR = q - p(:, k);
        dRU = p(:, k) - u;
        rBR = norm(dBR);
        rRU = norm(dRU);
        check_distance(rBR, sprintf('BS and RIS %d', k));
        check_distance(rRU, sprintf('RIS %d and UE', k));

        tau = (rBR + rRU)/c;
        nu = (-(pDot(:, k).')*dBR/rBR ...
            + (pDot(:, k) - uDot).'*dRU/rRU) / c;
        aSec(ell, k+1) = tau - (nu - eta)*T - delta;
    end
end
end


function check_distance(d, label)
if d < 1e-8
    error('Degenerate geometry: distance between %s is nearly zero.', label);
end
end


function Jprior = make_prior_fim(cfg, scale, idx)
Jprior = zeros(idx.nAll);
if ~cfg.risPrior.enabled
    return;
end
for j = idx.risPosition(:).'
    Jprior(j, j) = (scale(j)/cfg.risPrior.positionStdM)^2;
end
for j = idx.risVelocity(:).'
    Jprior(j, j) = (scale(j)/cfg.risPrior.velocityStdMps)^2;
end
end


function b = evaluate_at_snr(snrDb, signalPower, JdataUnitNoise, ...
        JpriorScaled, scale, idx, K)
sigma2 = signalPower / 10^(snrDb/10);
JgammaScaled = JdataUnitNoise/sigma2 + JpriorScaled;
JgammaScaled = symmetrize(JgammaScaled);

physical = 1:idx.nPhysical;
betaIdx = [idx.betaReal, idx.betaImag];
Jxx = JgammaScaled(physical, physical);
Jxb = JgammaScaled(physical, betaIdx);
Jbb = JgammaScaled(betaIdx, betaIdx);

% Equivalent FIM after eliminating the nuisance channel coefficients.
if rcond(Jbb) > 1e-12
    schurTerm = Jxb * (Jbb \ Jxb.');
else
    schurTerm = Jxb * pinv(Jbb) * Jxb.';
end
JeScaled = symmetrize(Jxx - schurTerm);

[CphysicalScaled, infoEquivalent] = stable_psd_inverse(JeScaled);
[CgammaScaled, infoFull] = stable_psd_inverse(JgammaScaled);

scalePhysical = scale(physical);
CphysicalSI = (scalePhysical * scalePhysical.') .* CphysicalScaled;
CgammaSI = (scale * scale.') .* CgammaScaled;
CphysicalSI = symmetrize(CphysicalSI);
CgammaSI = symmetrize(CgammaSI);

b = struct();
b.snrDb = snrDb;
b.noiseVariance = sigma2;
b.JgammaScaled = JgammaScaled;
b.JequivalentScaled = JeScaled;
b.CphysicalSI = CphysicalSI;
b.CgammaSI = CgammaSI;
b.infoEquivalent = infoEquivalent;
b.infoFull = infoFull;
b.stdPhysicalSI = sqrt(max(real(diag(CphysicalSI)), 0));
b.stdAllSI = sqrt(max(real(diag(CgammaSI)), 0));
b.uePEB_m = sqrt(max(real(trace(CphysicalSI(idx.uePosition, idx.uePosition))), 0));
b.ueVEB_mps = sqrt(max(real(trace(CphysicalSI(idx.ueVelocity, idx.ueVelocity))), 0));
b.etaVariance = real(CphysicalSI(idx.eta, idx.eta));
b.deltaVarianceSec2 = real(CphysicalSI(idx.delta, idx.delta));
b.etaStdPpm = sqrt(max(b.etaVariance, 0))*1e6;
b.deltaStdNs = sqrt(max(b.deltaVarianceSec2, 0))*1e9;
b.risPEB_m = zeros(K, 1);
b.risVEB_mps = zeros(K, 1);
for k = 1:K
    ip = idx.risPosition(:, k);
    iv = idx.risVelocity(:, k);
    b.risPEB_m(k) = sqrt(max(real(trace(CphysicalSI(ip, ip))), 0));
    b.risVEB_mps(k) = sqrt(max(real(trace(CphysicalSI(iv, iv))), 0));
end
b.betaRealStd = b.stdAllSI(idx.betaReal);
b.betaImagStd = b.stdAllSI(idx.betaImag);
end


function [C, info] = stable_psd_inverse(A)
A = symmetrize(A);
[V, D] = eig(A);
d = real(diag(D));
[d, order] = sort(d, 'ascend');
V = V(:, order);
dMax = max(abs(d));
tol = max(size(A)) * eps(max(dMax, 1));
positive = d > tol;

info.rank = sum(positive);
info.dimension = size(A, 1);
info.minEigenvalue = min(d);
info.maxEigenvalue = max(d);
if any(positive)
    info.conditionNumber = max(d(positive))/min(d(positive));
else
    info.conditionNumber = Inf;
end
info.usedPseudoInverse = ~all(positive);

if all(positive)
    C = V * diag(1./d) * V.';
else
    invD = zeros(size(d));
    invD(positive) = 1./d(positive);
    C = V * diag(invD) * V.';
    warning(['FIM rank is %d of %d. A Moore-Penrose inverse was used. ' ...
        'Some parameters are not locally identifiable under these settings.'], ...
        info.rank, info.dimension);
end
C = symmetrize(C);
end


function sweep = initialize_sweep(K, M, nSNR, nPhysical, nAll)
sweep.uePEB_m = zeros(1, nSNR);
sweep.ueVEB_mps = zeros(1, nSNR);
sweep.etaStdPpm = zeros(1, nSNR);
sweep.deltaStdNs = zeros(1, nSNR);
sweep.risPEB_m = zeros(K, nSNR);
sweep.risVEB_mps = zeros(K, nSNR);
sweep.betaRealStd = zeros(M, nSNR);
sweep.betaImagStd = zeros(M, nSNR);
sweep.stdPhysicalSI = zeros(nPhysical, nSNR);
sweep.stdAllSI = zeros(nAll, nSNR);
sweep.rankEquivalent = zeros(1, nSNR);
sweep.conditionEquivalent = zeros(1, nSNR);
end


function sweep = store_sweep_point(sweep, b, s)
sweep.uePEB_m(s) = b.uePEB_m;
sweep.ueVEB_mps(s) = b.ueVEB_mps;
sweep.etaStdPpm(s) = b.etaStdPpm;
sweep.deltaStdNs(s) = b.deltaStdNs;
sweep.risPEB_m(:, s) = b.risPEB_m;
sweep.risVEB_mps(:, s) = b.risVEB_mps;
sweep.betaRealStd(:, s) = b.betaRealStd;
sweep.betaImagStd(:, s) = b.betaImagStd;
sweep.stdPhysicalSI(:, s) = b.stdPhysicalSI;
sweep.stdAllSI(:, s) = b.stdAllSI;
sweep.rankEquivalent(s) = b.infoEquivalent.rank;
sweep.conditionEquivalent(s) = b.infoEquivalent.conditionNumber;
end


function T = make_parameter_table(b, idx, K)
labels = cell(idx.nAll, 1);
units = cell(idx.nAll, 1);
displayFactor = ones(idx.nAll, 1);
xyz = {'x', 'y', 'z'};

for d = 1:3
    labels{idx.uePosition(d)} = sprintf('UE position %s', xyz{d});
    units{idx.uePosition(d)} = 'm';
    labels{idx.ueVelocity(d)} = sprintf('UE velocity %s', xyz{d});
    units{idx.ueVelocity(d)} = 'm/s';
end
for k = 1:K
    for d = 1:3
        labels{idx.risPosition(d, k)} = sprintf('RIS %d position %s', k, xyz{d});
        units{idx.risPosition(d, k)} = 'm';
        labels{idx.risVelocity(d, k)} = sprintf('RIS %d velocity %s', k, xyz{d});
        units{idx.risVelocity(d, k)} = 'm/s';
    end
end
labels{idx.eta} = 'Relative frequency error eta';
units{idx.eta} = 'ppm';
displayFactor(idx.eta) = 1e6;
labels{idx.delta} = 'Clock offset delta_c';
units{idx.delta} = 'ns';
displayFactor(idx.delta) = 1e9;

for k = 1:(K+1)
    labels{idx.betaReal(k)} = sprintf('Re{beta_%d}', k-1);
    labels{idx.betaImag(k)} = sprintf('Im{beta_%d}', k-1);
    units{idx.betaReal(k)} = '1';
    units{idx.betaImag(k)} = '1';
end

varianceSI = max(real(diag(b.CgammaSI)), 0);
stdSI = sqrt(varianceSI);
stdReported = stdSI .* displayFactor;
T = table(labels, varianceSI, stdSI, stdReported, units, ...
    'VariableNames', {'Parameter', 'CRLBVarianceSI', 'StdLowerBoundSI', ...
    'StdLowerBoundReported', 'ReportedUnit'});
end


function T = make_sweep_table(sweep, K, M)
T = table(sweep.snrDb(:), sweep.uePEB_m(:), sweep.ueVEB_mps(:), ...
    sweep.etaStdPpm(:), sweep.deltaStdNs(:), ...
    'VariableNames', {'SNR_dB', 'UE_PEB_m', 'UE_VEB_mps', ...
    'eta_std_ppm', 'delta_std_ns'});
for k = 1:K
    T.(sprintf('RIS%d_PEB_m', k)) = sweep.risPEB_m(k, :).';
    T.(sprintf('RIS%d_VEB_mps', k)) = sweep.risVEB_mps(k, :).';
end
for k = 1:M
    T.(sprintf('beta%d_real_std', k-1)) = sweep.betaRealStd(k, :).';
    T.(sprintf('beta%d_imag_std', k-1)) = sweep.betaImagStd(k, :).';
end
end


function [fig1, fig2, fig3] = make_plots(cfg, b, sweep, idx, K, M)
snrDb = sweep.snrDb;

fig1 = figure('Color', 'w', 'Name', 'CRLB versus SNR');
subplot(2, 2, 1);
semilogy(snrDb, sweep.uePEB_m, 'o-', 'LineWidth', 1.5); hold on;
semilogy(snrDb, sweep.ueVEB_mps, 's-', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('Lower bound');
legend('UE PEB (m)', 'UE VEB (m/s)', 'Location', 'southwest');
title('UE bounds');

subplot(2, 2, 2);
for k = 1:K
    semilogy(snrDb, sweep.risPEB_m(k, :), 'o-', 'LineWidth', 1.3); hold on;
end
grid on; xlabel('SNR (dB)'); ylabel('RIS PEB (m)');
legend(make_ris_legend(K), 'Location', 'best'); title('RIS position bounds');

subplot(2, 2, 3);
for k = 1:K
    semilogy(snrDb, sweep.risVEB_mps(k, :), 's-', 'LineWidth', 1.3); hold on;
end
grid on; xlabel('SNR (dB)'); ylabel('RIS VEB (m/s)');
legend(make_ris_legend(K), 'Location', 'best'); title('RIS velocity bounds');

subplot(2, 2, 4);
semilogy(snrDb, sweep.etaStdPpm, 'd-', 'LineWidth', 1.5); hold on;
semilogy(snrDb, sweep.deltaStdNs, '^-', 'LineWidth', 1.5);
grid on; xlabel('SNR (dB)'); ylabel('Standard-deviation lower bound');
legend('eta (ppm)', 'delta_c (ns)', 'Location', 'southwest');
title('Synchronization bounds');

fig2 = figure('Color', 'w', 'Name', 'Component bounds at default SNR');
subplot(2, 2, 1);
bar(b.stdPhysicalSI(idx.uePosition)); grid on;
set(gca, 'XTickLabel', {'x', 'y', 'z'}); ylabel('Std. lower bound (m)');
title(sprintf('UE position at %.1f dB', cfg.snrDb));
subplot(2, 2, 2);
bar(b.stdPhysicalSI(idx.ueVelocity)); grid on;
set(gca, 'XTickLabel', {'v_x', 'v_y', 'v_z'}); ylabel('Std. lower bound (m/s)');
title('UE velocity components');
subplot(2, 2, 3);
bar(b.risPEB_m); grid on; ylabel('PEB (m)'); xlabel('RIS index');
title('Per-RIS position bounds');
subplot(2, 2, 4);
bar(b.risVEB_mps); grid on; ylabel('VEB (m/s)'); xlabel('RIS index');
title('Per-RIS velocity bounds');

fig3 = figure('Color', 'w', 'Name', 'Synchronization and channel bounds');
subplot(1, 2, 1);
bar([b.etaStdPpm, b.deltaStdNs]); set(gca, 'YScale', 'log'); grid on;
set(gca, 'XTickLabel', {'eta (ppm)', 'delta_c (ns)'});
ylabel('Reported standard-deviation lower bound');
title(sprintf('Synchronization at %.1f dB', cfg.snrDb));
subplot(1, 2, 2);
bar(0:M-1, [b.betaRealStd(:), b.betaImagStd(:)]); grid on;
xlabel('Path index k'); ylabel('Std. lower bound');
legend('Re{beta_k}', 'Im{beta_k}', 'Location', 'best');
title('Nuisance channel coefficients');
end


function labels = make_ris_legend(K)
labels = cell(1, K);
for k = 1:K
    labels{k} = sprintf('RIS %d', k);
end
end


function print_summary(b, cfg, K)
fprintf('\n============================================================\n');
fprintf('Hybrid CRLB summary at SNR = %.1f dB\n', cfg.snrDb);
fprintf('Geometry model: %s\n', cfg.geometryModel);
fprintf('------------------------------------------------------------\n');
fprintf('UE PEB                 : %.6g m\n', b.uePEB_m);
fprintf('UE VEB                 : %.6g m/s\n', b.ueVEB_mps);
for k = 1:K
    fprintf('RIS %d PEB / VEB        : %.6g m / %.6g m/s\n', ...
        k, b.risPEB_m(k), b.risVEB_mps(k));
end
fprintf('Std LB of eta          : %.6g ppm\n', b.etaStdPpm);
fprintf('Std LB of delta_c      : %.6g ns\n', b.deltaStdNs);
fprintf('Equivalent FIM rank    : %d / %d\n', ...
    b.infoEquivalent.rank, b.infoEquivalent.dimension);
fprintf('Equivalent FIM cond.   : %.6g (scaled coordinates)\n', ...
    b.infoEquivalent.conditionNumber);
fprintf('============================================================\n\n');
end


function save_outputs(results, parameterTable, sweepTable, fig1, fig2, fig3, cfg)
out = cfg.outputFolder;
if ~exist(out, 'dir')
    mkdir(out);
end
writetable(parameterTable, fullfile(out, 'default_parameter_bounds.csv'));
writetable(sweepTable, fullfile(out, 'snr_sweep_bounds.csv'));

resultsToSave = results;
resultsToSave = rmfield(resultsToSave, 'figureHandles');
save(fullfile(out, 'crlb_results.mat'), 'resultsToSave', '-v7');

if ~isempty(fig1)
    save_figure(fig1, fullfile(out, 'bounds_vs_snr.png'));
    save_figure(fig2, fullfile(out, 'component_bounds.png'));
    save_figure(fig3, fullfile(out, 'sync_and_channel_bounds.png'));
end
fprintf('Results saved in: %s\n', out);
end


function save_figure(fig, filename)
try
    exportgraphics(fig, filename, 'Resolution', 180);
catch
    saveas(fig, filename);
end
end


function A = symmetrize(A)
A = real((A + A.')/2);
end


function out = merge_struct(base, override)
% Recursive merge lets the user supply only the fields that should change.
out = base;
names = fieldnames(override);
for i = 1:numel(names)
    name = names{i};
    if isfield(out, name) && isstruct(out.(name)) && isstruct(override.(name))
        out.(name) = merge_struct(out.(name), override.(name));
    else
        out.(name) = override.(name);
    end
end
end
