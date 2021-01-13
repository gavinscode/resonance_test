clear; close all; clc

% Load data
nanosphere_reference

sizesToUse = [1 4 6];

diametersToCalc = (6:3:15)*10^-9;

% For 0th, 1st and 2nd modes
    % Set to -1 for no zeros (counted from 0)
zerosToCalc = [0 1 1];

zerosCols = ['r', 'b', 'g'];

figure;

% Plot linear function of diameter
subplot(1,2,1); hold on;
plot(nanocrystalSize_m*10^9, nanocrystalFreqResonance_hz/10^9, 'x')
plot(nanocrystalSize_m(sizesToUse)*10^9, nanocrystalFreqResonance_hz(sizesToUse)/10^9, 'o')

% Plot recipriocal function diameter
subplot(1,2,2); hold on;
plot(1./(nanocrystalSize_m*10^9), nanocrystalFreqResonance_hz/10^9, 'x')
plot(1./(nanocrystalSize_m(sizesToUse)*10^9), nanocrystalFreqResonance_hz(sizesToUse)/10^9, 'o')

%For each mode
for iMode = 0:2
    resonantFreqs = zeros(length(diametersToCalc), zerosToCalc(iMode+1)+1);
    
    if zerosToCalc(iMode+1) >= 0
        for jSize = 1:length(diametersToCalc)
            resonantFreqs(jSize,:) = calcualtesphereresonance(diametersToCalc(jSize)/2, ...
                'sph', iMode, zerosToCalc(iMode+1), CdSeVelocity_mps(1), CdSeVelocity_mps(2), 5*10^9, 10^6, 0);
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
ylim([100 300])
title('normal diameter')

subplot(1,2,2); hold on;
xlim([1/15 1/6])
title('reciprical diameter')
ylim([100 300])