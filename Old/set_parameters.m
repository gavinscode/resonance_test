LIGHT_SPEED = 299792458; % m/s

VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

% Freqs to cover - set to match range of q interpolation
testFrequncies_hz = (6:0.05:14)*10^9;

testFrequncies_rad = testFrequncies_hz*2*pi;

% From Yang 2016 - worked on H3N2
% Measured paramteres - for reference if changed in model
measuredResonance_hz = 8.22*10^9;

measuredResonance_rad = measuredResonance_hz*2*pi;

measuredBandwidth = 4.22*10^9;

measuredQualityFactor = 1.95;

measuredAbsorbtion = 0.21; % Max

% Correct reduced mass for 100 nm: 8.77*10^-21; ???
%%% Maybe this should be *9.8
reducedMass_kg = 14.5*1.6605*10^-21; % Convert from MDa

% Note this is density of virions in solution, not mass/density of particles
virionDensity = 7.5*10^14; % 1/m^3

channelLength_m = 1.25/1000;

% Correct radius for 161 MDa weight: 140 nm
diameter_m = 100*10^-9;

stressThreshold_pa = 0.141*10^6; %Pa

% Note replaced in charge distrubtion script
providedChargeDistribution = 1.16*10^7*1.602176634*10^-19;

% Model parameters
resonanceFrequency_hz = 8.22*10^9; 

resonantFrequency_rad = resonanceFrequency_hz*2*pi;

qualityFactor = 1.95;

systemSpring = resonantFrequency_rad^2*reducedMass_kg;

systemDamp = resonantFrequency_rad*reducedMass_kg/qualityFactor;

% Measured parameter from paper

% Absorbtion spectra
% From Fig. 3B of Yang 2016 
    % (y) absorbtion measured in imageJ for each Ghz
    % (x) frequncy measured for peak and zero
    % Peak is actually at 8.3??? Changed from 8.22
measuredFrequencies_hz = [6 7 8 8.3 9 10 11 12 13.25 14]*10^9;

measuredFrequencies_rad = measuredFrequencies_hz*2*pi;

measuredAbsorbtions = [9.7 16.4 21.1 21.5 20.1 13.8 6.4 1.7 0.1 0.5]/100;

% Inactivation results
% From Fig. 4b
inactivation1_fieldIntensity_vm = [273]; %v/m

inactivation1_result_ratio = [70.5 83.0 93.1 100 97.9 76.4 61.8 49.8];

inactivation1_frequncy_Hz = [6 7 8 8.4 9 10 11 12];

% From Fig. 5b
    % Point with intensity 273 inactivation 100% already included in previous
    % Mueasred to top and botom of points and averaged
inactivation2_fieldIntensity_vm = [68 87 171]; %v/m

inactivation2_result_ratio = [(4.9+6.9)/2 (36.7+39.0)/2 (62.8+65.2)/2];

inactivation2_frequncy_Hz = 8.4;

% from Ellison et al. 1996 pg 240 Pottel 1980, 25oC
    % Maybe 10% 50 oC change in temeprature on permitivity 
    % Use interpoalted eqns from p. 257
permitivityValues = [73.11, 73.19, 73.13, 72.71, 72.63, 72.25, 71.92, 71.50, ... 
    71.21, 70.34, 69.45, 69.74, 69.24, 68.78, 68.69, 67.86, 67.13, 66.77, ... 
    65.38, 63.76, 63.04, 61.74, 62.04, 61.10, 59.46, 58.49, 58.18, 55.78, ... 
    55.27, 54.03, 53.47, 52.63]; 

permitivityFrequenceis = [5.306, 5.323, 5.433, 5.536, 5.638, 5.853, 6.000, ...
    6.145, 6.300, 6.729, 6.850, 6.958, 7.267, 7.406, 7.681, 7.850, 8.243, ...
    8.579, 8.979, 9.516, 10.010, 10.150, 10.230, 10.450, 11.320, 11.730, ...
    12.000, 12.770, 13.140, 13.380, 13.820, 14.230];

