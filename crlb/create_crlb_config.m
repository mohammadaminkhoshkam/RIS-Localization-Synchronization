function cfg = create_crlb_config()
% Create the default CRLB simulation settings.

cfg.numSubcarriers = 129;
cfg.numSnapshots = 32;
cfg.subcarrierSpacingHz = 480e3;
cfg.snapshotIntervalSec = 5e-3;
cfg.centerObservationTimes = true;
cfg.ofdmEnergy = 1;
cfg.snrDb = 20;
cfg.snrSweepDb = -10:2.5:30;
cfg.randomSeed = 42;

% Node geometry
cfg.bsPositionM = [0; 0; 5];
cfg.uePositionM = [8; 7; 1.5];
cfg.ueVelocityMps = [1.2; -0.4; 0.1];
cfg.risPositionM = [0, 20, 18; ...
                    18, 0, 18; ...
                    3, 4, 6];
cfg.risVelocityMps = zeros(3, 3);

% Synchronization and channel
cfg.eta = 2e-6;
cfg.clockOffsetSec = 40e-9;
cfg.beta = [1.00*exp(1j*0.15); ...
            0.65*exp(-1j*0.80); ...
            0.50*exp(1j*1.40); ...
            0.40*exp(-1j*2.20)];

% RIS prior
cfg.risPrior.enabled = true;
cfg.risPrior.positionStdM = 0.10;
cfg.risPrior.velocityStdMps = 0.05;

cfg.geometryModel = 'reference-linearized';
cfg.makePlots = false;
cfg.verbose = false;
cfg.saveResults = false;
end
