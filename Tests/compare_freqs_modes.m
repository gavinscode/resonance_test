clear; close all; clc

% Load data
nanosphere_reference

sizesToUse = [1 4 6];

diametersToCalc = (6:2:16)*10^-9;

% For 0th, 1st and 2nd modes
    % Set to -1 for no zeros (counted from 0)
zerosToCalc = [0 1 1];

zerosCols = ['r', 'b', 'g'];

% Get core to shell mass fraction
coreFraction = zeros(length(sizesToUse),1);

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);
    
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    coreFraction(iSize) = coreMass/(coreMass + shellMass);
end

% Interp core fraction for test diameters, linear inside of range
coreFractionInterp = interp1(nanocrystalSize_m(sizesToUse), coreFraction, ...
    diametersToCalc, 'linear', -1);

% Extrapolate with nearest outside of range
toGet = find(coreFractionInterp == -1);

coreFractionInterp(toGet) = interp1(nanocrystalSize_m(sizesToUse), coreFraction, ...
    diametersToCalc(toGet), 'nearest', 'extrap');

figure;

% Plot linear function of diameter
subplot(1,2,1); hold on;
plot(nanocrystalSize_m*10^9, nanocrystalFreqResonance_hz/10^9, 'xk')
plot(nanocrystalSize_m(sizesToUse)*10^9, nanocrystalFreqResonance_hz(sizesToUse)/10^9, 'ok')

plot(nanocrystalSize_m*10^9, nanocrystal2ndFreqResonance_hz/10^9, 'xk')
plot(nanocrystalSize_m(sizesToUse)*10^9, nanocrystal2ndFreqResonance_hz(sizesToUse)/10^9, 'ok')

plot(nanocrystalSize_m*10^9, nanocrystal3rdFreqResonance_hz/10^9, 'xk')
plot(nanocrystalSize_m(sizesToUse)*10^9, nanocrystal3rdFreqResonance_hz(sizesToUse)/10^9, 'ok')

% Plot recipriocal function diameter
subplot(1,2,2); hold on;
plot(1./(nanocrystalSize_m*10^9), nanocrystalFreqResonance_hz/10^9, 'xk')
plot(1./(nanocrystalSize_m(sizesToUse)*10^9), nanocrystalFreqResonance_hz(sizesToUse)/10^9, 'ok')

plot(1./(nanocrystalSize_m*10^9), nanocrystal2ndFreqResonance_hz/10^9, 'xk')
plot(1./(nanocrystalSize_m(sizesToUse)*10^9), nanocrystal2ndFreqResonance_hz(sizesToUse)/10^9, 'ok')

plot(1./(nanocrystalSize_m*10^9), nanocrystal3rdFreqResonance_hz/10^9, 'xk')
plot(1./(nanocrystalSize_m(sizesToUse)*10^9), nanocrystal3rdFreqResonance_hz(sizesToUse)/10^9, 'ok')

%For each mode
for iMode = 0:2
    resonantFreqs = zeros(length(diametersToCalc), zerosToCalc(iMode+1)+1);
    
    if zerosToCalc(iMode+1) >= 0
        for jSize = 1:length(diametersToCalc)
            % Avg sound velocities based on core fractions
            avgLongVel = CdSeVelocity_mps(1)*coreFractionInterp(jSize) + ...
                CdTeVelocity_mps(1)*(1-coreFractionInterp(jSize));
            
            avgTransVel = CdSeVelocity_mps(2)*coreFractionInterp(jSize) + ...
                CdTeVelocity_mps(2)*(1-coreFractionInterp(jSize));
            
            resonantFreqs(jSize,:) = calcualtesphereresonance(diametersToCalc(jSize)/2, ...
                'sph', iMode, zerosToCalc(iMode+1), avgLongVel, avgTransVel, 5*10^9, 10^6, 0);
        end
        
        for jZero = 1:(zerosToCalc(iMode+1)+1)
            subplot(1,2,1); hold on;
            plot(diametersToCalc*10^9, resonantFreqs(:,jZero)/10^9, 'color', zerosCols(iMode+1))
            
            subplot(1,2,2); hold on;
            plot(1./(diametersToCalc*10^9), resonantFreqs(:,jZero)/10^9, 'color', zerosCols(iMode+1))
        end
    end
end

subplot(1,2,1); hold on;
xlim([6 15])
ylim([100 400])
title('normal diameter')

subplot(1,2,2); hold on;
xlim([1/15 1/6])
title('reciprical diameter')
ylim([100 400])