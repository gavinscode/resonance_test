close all;
clear;
clc

% Load and/or calculate values
set_parameters

set_charge_distribution

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

% Eqn 11
analyticalStress = 2*systemSpring*analyticalAmplitude/(0.58*pi*(diameter_m/2)^2);

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