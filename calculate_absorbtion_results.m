close all;
clear;
clc

% Load and/or calculate values
set_parameters

set_charge_distribution

%% Calculate analytical results
analyticalAmplitude_fixedCharge = zeros(length(testFrequncies_hz),1);

analyticalAbsorbtion_fixedCharge = zeros(length(testFrequncies_hz),1); 

analyticalAmplitude_varyingCharge = zeros(length(testFrequncies_hz),1);

analyticalAbsorbtion_varingCharge = zeros(length(testFrequncies_hz),1); 

% Note that eletric field amplitude is ignored for absorbtion as it cancels

for iFreq = 1:length(testFrequncies_hz)
    % Do first for fixed charge
    % Eqn 7
    analyticalAmplitude_fixedCharge(iFreq) = providedChargeDistribution./(reducedMass_kg*...
        sqrt((resonantFrequency_rad^2 - testFrequncies_rad(iFreq).^2).^2 + ...
        (resonantFrequency_rad*testFrequncies_rad(iFreq)/qualityFactor).^2));

    % Eqn 9
    analyticalPower = resonantFrequency_rad*...
        testFrequncies_rad(iFreq)^2 * reducedMass_kg * analyticalAmplitude_fixedCharge(iFreq)^2 / ...
        (2*qualityFactor);
    
    powerFlux = 0.5*sqrt(relativePermitivtyInterpolated(iFreq))*VACCUM_PERMITIVITY*LIGHT_SPEED;
    
    % Eqn 10
    theta_absorbtion = analyticalPower/powerFlux;
    
    % Eqn 13 - solution for absorbtion
    analyticalAbsorbtion_fixedCharge(iFreq) = 1-exp(-theta_absorbtion*virionDensity*channelLength_m);
    
    % Repeat for varying charge
    analyticalAmplitude_varyingCharge(iFreq) = qInterpolated(iFreq)./(reducedMass_kg*...
        sqrt((resonantFrequency_rad^2 - testFrequncies_rad(iFreq).^2).^2 + ...
        (resonantFrequency_rad*testFrequncies_rad(iFreq)/qualityFactor).^2));
    
    analyticalPower = resonantFrequency_rad*...
        testFrequncies_rad(iFreq)^2 * reducedMass_kg * analyticalAmplitude_varyingCharge(iFreq)^2 / ...
        (2*qualityFactor);
    
    powerFlux = 0.5*sqrt(relativePermitivtyInterpolated(iFreq))*VACCUM_PERMITIVITY*LIGHT_SPEED;
    
    theta_absorbtion = analyticalPower/powerFlux;
    
    analyticalAbsorbtion_varingCharge(iFreq) = 1-exp(-theta_absorbtion*virionDensity*channelLength_m);
end

% Eqn 11
analyticalStress_fixedCharge = 2*systemSpring*analyticalAmplitude_fixedCharge/(0.58*pi*(diameter_m/2)^2);

analyticalStress_varyingCharge = 2*systemSpring*analyticalAmplitude_varyingCharge/(0.58*pi*(diameter_m/2)^2);

%% Plot results

figure; 
yyaxis left; hold on
plot(testFrequncies_hz/10^9, analyticalStress_varyingCharge/1000, 'k')

ylim([0 2])

ylabel('Stress (kPa)'); 

set(gca,'ycolor','k')

set(gca, 'YTick', 0:0.5:2);

yyaxis right; hold on
plot(testFrequncies_hz/10^9, analyticalAbsorbtion_varingCharge*100, 'b')

ylabel('Absorbtion (%)')

xlabel('Frequency (GHz)')

set(gca,'TickDir','out');

set(gcf, 'Position', [50 50 500 500/1.61]);
set(gca,'ycolor', 'b')

figure; hold on
plot(testFrequncies_hz/10^9, analyticalAbsorbtion_fixedCharge*100, 'r')

plot(testFrequncies_hz/10^9, analyticalAbsorbtion_varingCharge*100, 'b')

legend('Fixed Charge', 'Varying charge')

ylabel('Absorbtion (%)')

xlabel('Frequency (GHz)')

set(gca,'TickDir','out');

set(gcf, 'Position', [50 50 500 500/1.61]);

%% Get peak 
% [peakAbsorbtion, peakInd] = max(analyticalAbsorbtion);
% 
% plot(testFrequncies_hz(peakInd)/10^9, analyticalAbsorbtion(peakInd), 'o')
% 
% % Indicate expected resonance
% [~, resonanceInd] = min(abs(testFrequncies_hz - measuredResonance_hz));
% 
% plot(testFrequncies_hz(resonanceInd)/10^9, analyticalAbsorbtion(resonanceInd), 'x');
% 
% % Indicate bandwidth
% upperBandInd = find(analyticalAbsorbtion(peakInd:end) < measuredAbsorbtion/2);
% 
% upperBandInd = upperBandInd(1) + peakInd - 1;
% 
% lowerBandInd = find(analyticalAbsorbtion(1:peakInd) < measuredAbsorbtion/2);
% 
% if ~isempty(lowerBandInd)
%    lowerBandInd = lowerBandInd(end);
%    
%    plot(testFrequncies_hz([lowerBandInd upperBandInd])/10^9, analyticalAbsorbtion([lowerBandInd upperBandInd]), '-')
%    
%    absorbtionBandwdith = (testFrequncies_hz(upperBandInd) - testFrequncies_hz(lowerBandInd))
%    
%    absorbtionQualityFactor = testFrequncies_hz(peakInd)/absorbtionBandwdith
% end
% 
% subplot(1,3,3); hold on
% plot(testFrequncies_hz([1 end])/10^9, [1 1]*stressThreshold_pa)
% 
% for iField = 1:length(fieldIntensity)
%     plot(testFrequncies_hz/10^9, analyticalStress*fieldIntensity(iField))
% end