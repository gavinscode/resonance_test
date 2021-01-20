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

    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    
    QEstimate(iSize) = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);
    
    
    resonance_rad = nanocrystalFreqResonance_hz(sizeIndex)*2*pi;

    systemSpring = resonance_rad^2*reducedMass;

    systemDamp = resonance_rad*reducedMass/QEstimate(iSize);
    
    qEstimate(iSize) = sqrt(nanocrystalThetaEx_m2(sizeIndex) * resonance_rad * ...
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
ylim([0 500]); 
xlim([0 15]);

% Note this is e as function of nanometers...

% Use fit 2nd order polynomial
opts = fitoptions('poly2', 'Lower', [-Inf 0 0], 'Upper', [Inf 0 0]);
fittedQuad = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'poly2', opts)

quadH = plot(fittedQuad,'b');

opts = fitoptions('exp1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedExp = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'exp1', opts)

expH = plot(fittedExp, 'g');

opts = fitoptions('power1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedPower = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'power1', opts)

powerH = plot(fittedPower, 'r');

legend([quadH, expH, powerH], {'Quad', 'Exp', 'Power'})