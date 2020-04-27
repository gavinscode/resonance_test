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

%To compare to absorbtion CS
[sqrt(theta_abs*resonantFrequency_rad*reducedMass_kg*LIGHT_SPEED*...
    sqrt(relativePermitivity)*vaccumPermitivity/qualityFactor)  chargeDistribution]

% Eqn 12 - works
radius_m = 50*10^-9;

PThreshold_pa = 0.141*10^6;

testFrequency_hz = [6 8 10]*10^9;

testFrequency_rad = testFrequency_hz*2*pi;

energyThreshold = PThreshold_pa*pi*radius_m^2*sqrt(...
    reducedMass_kg^2*(resonantFrequency_rad^2-testFrequency_rad.^2).^2 + ...
    ((resonantFrequency_rad*reducedMass_kg/qualityFactor)^2)*testFrequency_rad.^2) / ...
    (3.45*chargeDistribution*reducedMass_kg*resonantFrequency_rad^2)

% Need to integrate over one full cycle to get Pabsorbtion
%%% Compare numerical to theoretical
% Eqn 7
driveEnergy = 1;
displacmentAmplitued = chargeDistribution*driveEnergy./(reducedMass_kg*...
    sqrt((resonantFrequency_rad^2-testFrequency_rad.^2).^2 + ...
    (resonantFrequency_rad.*testFrequency_rad/qualityFactor).^2))

% Eqn 9
integratedPower = resonantFrequency_rad.*testFrequency_rad.^2*reducedMass_kg.*...
    displacmentAmplitued.^2/(2*qualityFactor)

% Numerical simulation (Euler method - simple) with forcing function
% For sine driving signal
sineOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*targetSignal(i-1)-systemSpring*sineOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    sineOutput(i) = sineOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

max(sineOutput)/displacmentAmplitued(2)

% For pulsed driving signal
pulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*pulsedSignal(i-1)-systemSpring*pulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    pulsedOutput(i) = pulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For pulsed lightwave driving signal
lightPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*lightPulseSignal(i-1)-systemSpring*lightPulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    lightPulsedOutput(i) = lightPulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

max(lightPulsedOutput)/displacmentAmplitued(2)

% For beating driving signal
beatingOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*beatingSignal(i-1)-systemSpring*beatingOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    beatingOutput(i) = beatingOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For femtopulsed driving signal
femtoPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*femtoPulsedSignal(i-1)-systemSpring*femtoPulsedOutput(i-1)-...
        systemDamp*dY)/reducedMass_kg;
    
    femtoPulsedOutput(i) = femtoPulsedOutput(i-1) + dY*samplingPeriod_s;
    
    dY = dY + ddY*samplingPeriod_s;
end

% For lightwave femtopulsed driving signal
lighFemtoPulsedOutput = zeros(sampleLength,1);

dY = 0;

for i = 2:sampleLength
    
    ddY = (chargeDistribution*lightFemtoPulsedSignal(i-1)-systemSpring*lighFemtoPulsedOutput(i-1)-...
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