%%% Should confirm with more accurate ODE solver 

% Influenza parameters from Yang 2016
resonanceFrequency_hz = 8.22*10^9; 

resonantFrequency_rad = resonanceFrequency_hz*2*pi;

qualityFactor = 1.92;

reducedMass_kg = 14.5*1.6605*10^-21; % Convert from MDa to kg

systemSpring = resonantFrequency_rad^2*reducedMass_kg;

systemDamp = resonantFrequency_rad*reducedMass_kg/qualityFactor;

% Taking E to denote elementry charge in C
chargeDistribution = 1.16*10^7*1.602176634*10^-19;

% Eqn 14 - works
relativePermitivity = 67.13; %At 8.2 GHz

vaccumPermitivity = 8.854187817*10^-12; %C^2/(N.M^2)

theta_abs = 2.5*10^-13; %M^2

%To compare to charge distribution
sqrt(theta_abs*resonantFrequency_rad*reducedMass_kg*LIGHT_SPEED*...
    sqrt(relativePermitivity)*vaccumPermitivity/qualityFactor) / chargeDistribution;

% Eqn 12 - works - only correct for sine forcing function
radius_m = 50*10^-9;

PThreshold_pa = 0.141*10^6;

testFrequency_hz = [6 8.22 10]*10^9;

testFrequency_rad = testFrequency_hz*2*pi;

energyThreshold = PThreshold_pa*pi*radius_m^2*sqrt(...
    reducedMass_kg^2*(resonantFrequency_rad^2-testFrequency_rad.^2).^2 + ...
    ((resonantFrequency_rad*reducedMass_kg/qualityFactor)^2)*testFrequency_rad.^2) / ...
    (3.45*chargeDistribution*reducedMass_kg*resonantFrequency_rad^2)

% Need to integrate over one full cycle to get Pabsorbtion
%%% Compare numerical to theoretical
% Eqn 7
driveEnergy = energyThreshold(2);
displacmentAmplitued = chargeDistribution*driveEnergy./(reducedMass_kg*...
    sqrt((resonantFrequency_rad^2-testFrequency_rad(2).^2).^2 + ...
    (resonantFrequency_rad.*testFrequency_rad(2)/qualityFactor).^2));

% Eqn 9
integratedPower = resonantFrequency_rad.*testFrequency_rad(2).^2*reducedMass_kg.*...
    displacmentAmplitued.^2/(2*qualityFactor);

% Approximate for carrier frequncy
% Scale up absorbtion based on absorbtion in water at carrier and resonant frequncy
resonance_absorbtion = (5.5536 + 5.2881)/2;

% Choose correct one for frequncy
%carrier_absorbtion = 1.071; % For 1259 nm - gain here 30.9 - stress ratio  0.000483

%carrier_absorbtion = 5.568; % For 1671 nm - gain here 61.1 - stress ratio 0.0017

carrier_absorbtion = 19.278; % For 2208 nm - gain here 98.9 - stress ratio 0.0027

theta_abs_scaled = -1/(1.25/1000 * 7.5*10^14)*log((1-0.21)^ ...
    (carrier_absorbtion/resonance_absorbtion));

carrierFrequency_rad = carrierFrequency_hz*2*pi;

% 1.77 is relative permittivity of light in water
chargeDistributionTemp = sqrt(theta_abs_scaled*carrierFrequency_rad*reducedMass_kg*LIGHT_SPEED*...
    sqrt(1.77)*vaccumPermitivity/qualityFactor);

chargeDistributionTemp/chargeDistribution

warning('Modified charge distrubtion used for all carrier waves')

% Numerical simulation (Euler method - simple) with forcing function
% For sine driving signal
sineOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*driveEnergy*targetSignal(i-1)-systemSpring*sineOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    sineOutput(i) = sineOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

sine_amp = (max(sineOutput)-min(sineOutput))/2;
sine_amp/displacmentAmplitued

% For pulsed driving signal
pulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*driveEnergy*pulsedSignal(i-1)-systemSpring*pulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    pulsedOutput(i) = pulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For pulsed lightwave driving signal
lightPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistributionTemp*driveEnergy*lightPulseSignal(i-1)-systemSpring*lightPulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    lightPulsedOutput(i) = lightPulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

amp_pulsed = (max(lightPulsedOutput)-min(lightPulsedOutput))/2;

amp_pulsed/displacmentAmplitued

%%% Test this against stress on normal...
stress_pulsed = 2*systemSpring*amp_pulsed/(0.58*pi*radius_m^2);

stress_pulsed/PThreshold_pa

% For beating driving signal
beatingOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistributionTemp*driveEnergy*beatingSignal(i-1)-systemSpring*beatingOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    beatingOutput(i) = beatingOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For femtopulsed driving signal
femtoPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistributionTemp*driveEnergy*femtoPulsedSignal(i-1)-systemSpring*femtoPulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    femtoPulsedOutput(i) = femtoPulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For lightwave femtopulsed driving signal
lighFemtoPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*driveEnergy*lightFemtoPulsedSignal(i-1)-systemSpring*lighFemtoPulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    lighFemtoPulsedOutput(i) = lighFemtoPulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end



% Plotting
figure; subplot(2,3,1); hold on 

plot(time_s, sineOutput*10^12); title('Sine drive');

subplot(2,3,2); hold on 

plot(time_s, pulsedOutput*10^12); title('Pulse drive')

subplot(2,3,3); hold on 

plot(time_s, femtoPulsedOutput*10^12); title('Femto-pulse drive')

subplot(2,3,4); hold on

plot(time_s, lightPulsedOutput*10^12); title('Pulsed lightwave drive')

subplot(2,3,5); hold on

plot(time_s, beatingOutput*10^12); title('Beating drive')

subplot(2,3,6); hold on

plot(time_s, lighFemtoPulsedOutput*10^12); title('Femtopulsed lightwave drive')