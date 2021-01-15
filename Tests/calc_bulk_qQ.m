clear; close all; clc

% Load data
nanosphere_reference

% Interpolate Core Size
sizesToUse = [1 4 6];

coresToGet = find(isnan(nanocrystalCore_m));

nanocrystalCore_m(coresToGet) = interp1(nanocrystalSize_m(sizesToUse), nanocrystalCore_m(sizesToUse), ...
    nanocrystalSize_m(coresToGet), 'linear');

% Calc q and Q values from original to get rough idea on relationship
sizesToUse = 1:6;    

qEstimate = zeros(length(sizesToUse),1);

QEstimate = zeros(length(sizesToUse),1);

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    QEstimate(iSize) = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);

    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    qEstimate(iSize) = sqrt(nanocrystalThetaEx_m2(sizeIndex) * nanocrystalFreqResonance_hz(sizeIndex) * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED./QEstimate(iSize)); 
end

figure;
subplot(1,2,1)
plot(nanocrystalSize_m(sizesToUse)*10^9, QEstimate, 'x-')
title('Quality factor');
ylim([0 15]); xlim([5 15]);

subplot(1,2,2); hold on
plot(nanocrystalSize_m(sizesToUse)*10^9, qEstimate/(1.602176634*10^-19), 'x-')
title('Charge (in e)');
ylim([0 150]); 
xlim([0 15]);

% Use fit 2nd order polynomial
opts = fitoptions('poly2', 'Lower', [-Inf 0 0], 'Upper', [Inf Inf 0]);
fittedQuad = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'poly2', opts)

quadH = plot(fittedQuad,'b');

opts = fitoptions('exp1', 'Lower', [1 -Inf], 'Upper', [1 Inf]);
fittedExp = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'exp1', opts)

expH = plot(fittedExp, 'g');

opts = fitoptions('power1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedPower = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'power1', opts)

powerH = plot(fittedPower, 'r');

legend([quadH, expH, powerH], {'Quad', 'Exp', 'Power'})