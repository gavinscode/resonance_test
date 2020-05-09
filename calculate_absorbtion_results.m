% From  Yang 2016, Figure 3B 
close all;
clear;
clc

LIGHT_SPEED = 299792458; % m/s

VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

% Freqs to cover - set to match range of q interpolation
testFrequncies_hz = (6:0.05:13)*10^9;

testFrequncies_rad = testFrequncies_hz*2*pi;

% Measured paramteres - for reference if changed in model
measuredResonance_hz = 8.22*10^9;

measuredResonance_rad = measuredResonance_hz*2*pi;

measuredBandwidth = 4.22*10^9;

measuredQualityFactor = 1.95;

measuredAbsorbtion = 0.21; % Max

reducedMass_kg = 14.5*1.6605*10^-21; % Convert from MDa to kg

virionDensity = 7.5*10^14; % 1/m^3

channelLength_m = 1.25/1000;

radius_m = 50*10^-9;

stressThreshold_pa = 0.141*10^6;

assumedChargeDistribution = 1.16*10^7*1.602176634*10^-19; % At 8.2 Ghz, 

% Model parameters
resonanceFrequency_hz = 8.22*10^9; 

resonantFrequency_rad = resonanceFrequency_hz*2*pi;

qualityFactor = 1.95;

systemSpring = resonantFrequency_rad^2*reducedMass_kg;

systemDamp = resonantFrequency_rad*reducedMass_kg/qualityFactor;

fieldIntensity = [68 87 171 274]; %v/m

%% To determine q (charge distribution) values
calcFrequencies_hz = [6 8.22 10.5 11.25 13]*10^9;

calcFrequencies_rad = calcFrequencies_hz*2*pi;

calcAbsorbtions = [0.1 0.21 0.1 0.05 0.001];

% from Ellison et al. 1996 pg 240 Pottel 1980, 25oC
% Exact freqs used are: 
% 6.000 8.243 10.450 11.320 13.140
calcRelativePermitivity = [71.92 67.13 61.10 59.46 55.27]; % Relative, 

% Calcualte values for q
% Eqn 7: q omitted as it is solved for and E is ommited as it will cancel later
calcAmplitudeWOq = 1./(reducedMass_kg*...
        sqrt((measuredResonance_rad^2 - calcFrequencies_rad.^2).^2 + ...
        (measuredResonance_rad*calcFrequencies_rad/measuredQualityFactor).^2));
 
% Eqn 9 - now omitting q^2   
calcAvgAbsorbtion = measuredResonance_rad.*calcFrequencies_rad.^2*...
    reducedMass_kg.*calcAmplitudeWOq.^2/(2*measuredQualityFactor);

% E^2 is ommited as will cancel later
calcPowerFlux = 0.5*sqrt(calcRelativePermitivity)*VACCUM_PERMITIVITY*LIGHT_SPEED;

% Eqn 10, q^2 is omitted as it will be solved for
calcThetaAbs1 = calcAvgAbsorbtion./calcPowerFlux;

% Eqn 13
calcThetaAbs2 = -1/(virionDensity*channelLength_m).*log(1-calcAbsorbtions);

% Calculate Q from equvielent between Theta Abs
calcQ = sqrt(calcThetaAbs2./calcThetaAbs1);

% Relative to orignal charge distrubiton
calcQ/assumedChargeDistribution;

% Compare to eqn 14 - simplified by assuming freq at resonance
sqrt(calcThetaAbs2(2)*measuredResonance_rad*reducedMass_kg*...
    sqrt(calcRelativePermitivity(2))*VACCUM_PERMITIVITY*LIGHT_SPEED...
    /measuredQualityFactor)/assumedChargeDistribution;

% Interpolate for charge distribution up to top freq calculated
qInterpolated = interp1(calcFrequencies_hz, calcQ, testFrequncies_hz, ...
    'linear','extrap');

% Quick plot to test results
figure; hold on; 

plot(testFrequncies_hz/10^9, qInterpolated/assumedChargeDistribution)

plot(calcFrequencies_hz/10^9, calcQ/assumedChargeDistribution, 'ro')

% Interpolate permitivity
relativePermitivtyInterpolated = interp1(calcFrequencies_hz, calcRelativePermitivity, ...
    testFrequncies_hz, 'linear','extrap');
%% Messy test for IR - but x17,000 gain!!! - probably not realistic
%%% May need to use resonance frequncy corresponding to higher mode, then 
% other amplitudes will be larger.

IR_WL = 1671*10^-9; % Middle 3rd water window
IR_freq_Rad = LIGHT_SPEED/IR_WL*2*pi;
% Assume absorbtion is simliar, so ThetaAbs2 will match as well
IR_permit = 1.77;

% This ends up super small - seems fair
IR_Amp = 1./(reducedMass_kg*...
        sqrt((measuredResonance_rad^2 - IR_freq_Rad.^2).^2 + ...
        (measuredResonance_rad*IR_freq_Rad/measuredQualityFactor).^2));

% This ends up relatively small - seems fair
IR_AvgAbs = measuredResonance_rad.*IR_freq_Rad.^2*...
    reducedMass_kg.*IR_Amp.^2/(2*measuredQualityFactor);

% This ends up a bit smaller - expected due to permitivity
IR_Flux = 0.5*sqrt(IR_permit)*VACCUM_PERMITIVITY*LIGHT_SPEED;

% This ends up relatively small - seems fair
IR_ThetaAbs1 = IR_AvgAbs/IR_Flux;

% Segelstein gives similar absorbtion coefficient for water between third water window and 8 Ghz, 
% Does this mean absorbtion ratio is similiar? If so, is cross-section equal?
IR_ThetaAbs2 = 2.5*10^-13;

% Because thetaAbs1 decrease but thetaAbs2 stays constant, q increases a lot...
IR_calcQ = sqrt(IR_ThetaAbs2./IR_ThetaAbs1);
    
IR_calcQ/assumedChargeDistribution;

%% Calculate analytical results
analyticalAmplitude = zeros(length(testFrequncies_hz),1);

analyticalAbsorbtion = zeros(length(testFrequncies_hz),1); 

% Note that eletric field amplitude is ignored for absorbtion as it cancels

for iFreq = 1:length(testFrequncies_hz)
    % Eqn 7
    analyticalAmplitude(iFreq) = qInterpolated(iFreq)./(reducedMass_kg*...
        sqrt((resonantFrequency_rad^2 - testFrequncies_rad(iFreq).^2).^2 + ...
        (resonantFrequency_rad*testFrequncies_rad(iFreq)/qualityFactor).^2));

    % Eqn 9
    analyticalPower = resonantFrequency_rad*...
        testFrequncies_rad(iFreq)^2 * reducedMass_kg * analyticalAmplitude(iFreq)^2 / ...
        (2*qualityFactor);
    
    powerFlux = 0.5*sqrt(relativePermitivtyInterpolated(iFreq))*VACCUM_PERMITIVITY*LIGHT_SPEED;
    
    % Eqn 10
    theta_absorbtion = analyticalPower/powerFlux;
    
    % Eqn 13 - solution for absorbtion
    analyticalAbsorbtion(iFreq) = 1-exp(-theta_absorbtion*virionDensity*channelLength_m);
end

analyticalStress = 2*systemSpring*analyticalAmplitude/(0.58*pi*radius_m^2);

%% Plot results

figure;

subplot(1,3,1); 
plot(testFrequncies_hz/10^9, analyticalAmplitude*10^12)

subplot(1,3,2); hold on
plot(testFrequncies_hz/10^9, analyticalAbsorbtion)

% Get peak 
[peakAbsorbtion, peakInd] = max(analyticalAbsorbtion);

plot(testFrequncies_hz(peakInd)/10^9, analyticalAbsorbtion(peakInd), 'o')

% Indicate expected resonance
[~, resonanceInd] = min(abs(testFrequncies_hz - measuredResonance_hz));

plot(testFrequncies_hz(resonanceInd)/10^9, analyticalAbsorbtion(resonanceInd), 'x');

% Indicate bandwidth
upperBandInd = find(analyticalAbsorbtion(peakInd:end) < measuredAbsorbtion/2);

upperBandInd = upperBandInd(1) + peakInd - 1;

lowerBandInd = find(analyticalAbsorbtion(1:peakInd) < measuredAbsorbtion/2);

if ~isempty(lowerBandInd)
   lowerBandInd = lowerBandInd(end);
   
   plot(testFrequncies_hz([lowerBandInd upperBandInd])/10^9, analyticalAbsorbtion([lowerBandInd upperBandInd]), '-')
   
   absorbtionBandwdith = (testFrequncies_hz(upperBandInd) - testFrequncies_hz(lowerBandInd))
   
   absorbtionQualityFactor = testFrequncies_hz(peakInd)/absorbtionBandwdith
end

subplot(1,3,3); hold on
plot(testFrequncies_hz([1 end])/10^9, [1 1]*stressThreshold_pa)

for iField = 1:length(fieldIntensity)
    plot(testFrequncies_hz/10^9, analyticalStress*fieldIntensity(iField))
end