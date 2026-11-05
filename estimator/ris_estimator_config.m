function cfg = ris_estimator_config()
% Default settings for the hybrid MAP/ML estimator.

cfg.c = 299792458;
cfg.N = 129;
cfg.L = 32;
cfg.deltaF = 480e3;
cfg.snapshotInterval = 5e-3;
cfg.centerTimes = true;
cfg.energy = 1;
cfg.snrDb = 20;
cfg.randomSeed = 42;
cfg.geometryModel = 'reference-linearized';

cfg.q = [3; 4; 5];

cfg.true.u0 = [8; 10; 1.5];
cfg.true.uDot = [1.2; -0.4; 0.1];
cfg.true.p0 = [0, 20, 18; ...
               18, 0, 18; ...
               3, 4, 6];
cfg.true.pDot = zeros(3, 3);
cfg.true.eta = 2e-6;
cfg.true.deltaC = 40e-9;
cfg.true.beta = [1.00*exp(1j*0.15); ...
                 0.65*exp(-1j*0.80); ...
                 0.50*exp(1j*1.40); ...
                 0.40*exp(-1j*2.20)];

priorOffset = [0.04, -0.03, 0.05; ...
              -0.02, 0.05, -0.04; ...
               0.03, -0.02, 0.02];
cfg.prior.enabled = true;
cfg.prior.pMean = cfg.true.p0 + priorOffset;
cfg.prior.pDotMean = zeros(size(cfg.true.pDot));
cfg.prior.positionStd = 0.10;
cfg.prior.velocityStd = 0.05;

cfg.initial.u0 = [7.8; 6.8; 1.7];
cfg.initial.uDot = [1.0; -0.3; 0.0];
cfg.initial.p0 = cfg.prior.pMean;
cfg.initial.pDot = cfg.prior.pDotMean;
cfg.initial.eta = 1e-6;
cfg.initial.deltaC = 30e-9;

cfg.bounds.u0Lower = [0; 0; 0];
cfg.bounds.u0Upper = [25; 25; 10];
cfg.bounds.uDotLower = [-5; -5; -5];
cfg.bounds.uDotUpper = [5; 5; 5];
cfg.bounds.p0Lower = cfg.prior.pMean - 0.5;
cfg.bounds.p0Upper = cfg.prior.pMean + 0.5;
cfg.bounds.pDotLower = -0.3*ones(size(cfg.true.pDot));
cfg.bounds.pDotUpper = 0.3*ones(size(cfg.true.pDot));
cfg.bounds.etaLower = -20e-6;
cfg.bounds.etaUpper = 20e-6;
cfg.bounds.deltaLower = -300e-9;
cfg.bounds.deltaUpper = 300e-9;

cfg.finiteDifference.positionStep = 1e-3;
cfg.finiteDifference.velocityStep = 1e-2;

cfg.optimizer.numStarts = 8;
cfg.optimizer.candidateMultiplier = 4;
cfg.optimizer.maxIterations = 600;
cfg.optimizer.maxFunctionEvaluations = 30000;
cfg.optimizer.optimalityTolerance = 1e-7;
cfg.optimizer.stepTolerance = 1e-8;
cfg.optimizer.display = 'final';
cfg.optimizer.startStd.u0 = 0.25;
cfg.optimizer.startStd.uDot = 0.15;
cfg.optimizer.startStd.p0 = 0.03;
cfg.optimizer.startStd.pDot = 0.015;
cfg.optimizer.startStd.eta = 0.15e-6;
cfg.optimizer.startStd.deltaC = 2e-9;

cfg.optimizer.syncSearch.enabled = true;
cfg.optimizer.syncSearch.etaHalfSpan = 4e-6;
cfg.optimizer.syncSearch.deltaHalfSpan = 80e-9;
cfg.optimizer.syncSearch.coarseEtaPoints = 25;
cfg.optimizer.syncSearch.coarseDeltaPoints = 33;
cfg.optimizer.syncSearch.fineEtaPoints = 13;
cfg.optimizer.syncSearch.fineDeltaPoints = 13;

cfg.optimizer.localBounds.enabled = true;
cfg.optimizer.localBounds.u0Radius = [3; 3; 3];
cfg.optimizer.localBounds.uDotRadius = [2; 2; 2];
cfg.optimizer.localBounds.p0Radius = 0.35;
cfg.optimizer.localBounds.pDotRadius = 0.20;

cfg.optimizer.maxNormalizedDataCost = 2.5;

cfg.saveResults = true;
cfg.outputFolder = 'hybrid_map_results';
end
