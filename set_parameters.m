LIGHT_SPEED = 299792458; % m/s

VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

% Freqs to cover - set to match range of q interpolation
testFrequncies_hz = (6:0.05:13)*10^9;

testFrequncies_rad = testFrequncies_hz*2*pi;

% From Yang 2016
% Measured paramteres - for reference if changed in model
measuredResonance_hz = 8.22*10^9;

measuredResonance_rad = measuredResonance_hz*2*pi;

measuredBandwidth = 4.22*10^9;

measuredQualityFactor = 1.95;

measuredAbsorbtion = 0.21; % Max

% Correct reduced mass for 100 nm: 8.77*10^-21;
reducedMass_kg = 14.5*1.6605*10^-21; % Convert from MDa

% Note this is density of virions in solution, not mass/density of particles
virionDensity = 7.5*10^14; % 1/m^3

channelLength_m = 1.25/1000;

% Correct radius for 161 MDa weight: 140 nm
diameter_m = 100*10^-9;

stressThreshold_pa = 0.141*10^6;

providedChargeDistribution = 1.16*10^7*1.602176634*10^-19; % At 8.2 Ghz, 

fieldIntensity = [68 87 171 274]; %v/m

% Model parameters
resonanceFrequency_hz = 8.22*10^9; 

resonantFrequency_rad = resonanceFrequency_hz*2*pi;

qualityFactor = 1.95;

systemSpring = resonantFrequency_rad^2*reducedMass_kg;

systemDamp = resonantFrequency_rad*reducedMass_kg/qualityFactor;