clear; close all; clc

% Load data
nanosphere_reference

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

% Calc 1 to 3rd order polynomials and compare

cols = ['r', 'b', 'g'];

for iOrder = 1:3
    qPoly = polyfit([nanocrystalSize_m], [qEstimate'], iOrder);

    qValues = polyval(qPoly, (1:15)/10^9);
    
    plot((1:15), qValues/(1.602176634*10^-19), ':', 'color', cols(iOrder))
    
end

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    
    
% Unkowns are q and Q - firstly determine just using abs. and bandwidth
frequencyRange = (100:300)*10^9*2*pi;

% Blurring filter from source bandwidth

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    coreFraction = coreMass/(coreMass + shellMass);
    
    avgLongVel = CdSeVelocity_mps(1)*coreFraction + ...
                CdTeVelocity_mps(1)*(1-coreFraction);
            
    avgTransVel = CdSeVelocity_mps(2)*coreFraction + ...
        CdTeVelocity_mps(2)*(1-coreFraction);
    
     % Get calcualted parameter values
        % Note, all frequenciesa are in radians
    resonance = calcualtesphereresonance(nanocrystalSize_m(sizeIndex)/2, ...
            'sph', 1, 0, avgLongVel, avgTransVel, 5*10^9, 10^6, 1)*2*pi;
    
    %%% Q is an unknown, unsure how it should vary with size
    systemQ = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);
    
    systemSpring = resonance^2*reducedMass;
    
    systemDamp = resonance*reducedMass/systemQ;

    %%% q is an unknown and should vary with the square of size 
    qToUse = sqrt(nanocrystalThetaEx_m2(sizeIndex) * nanocrystalFreqResonance_hz(sizeIndex) * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED/systemQ);  
    
    analyticalAbsorbtion = zeros(length(frequencyRange), 1);
    
    extinctionCrossSection = zeros(length(frequencyRange), 1);
    
    if freqResolution_Ghz(sizeIndex) == 50
        % convert freq resolution (when FWHM) to sigma
        sigma = freqResolution_Ghz(sizeIndex)/(2*sqrt(2*log(2)));
        
        sourceSpectra = 1/(sigma*sqrt(2*pi))*exp(-((1:200)-100).^2/...
            (2*(sigma)^2));
    else
        error('No conversion for resolution to sigma')
    end
    
    for jFreq = 1:length(frequencyRange)
        % Eqn 7
        analyticalAmplitude = qToUse./(reducedMass*...
            sqrt((resonance^2 - frequencyRange(jFreq)^2)^2 + ...
            (resonance*frequencyRange(jFreq)/systemQ).^2));

        % Eqn 9
        analyticalPower = resonance * ...
            frequencyRange(jFreq)^2 * reducedMass * ...
            analyticalAmplitude^2 / (2*systemQ);

        powerFlux = 0.5*VACCUM_PERMITIVITY*LIGHT_SPEED;

        % Eqn 10
        extinctionCrossSection(jFreq) = analyticalPower/powerFlux;
        
        % Eqn 13 - solution for absorbtion
        analyticalAbsorbtion(jFreq) = (1-exp(-extinctionCrossSection(jFreq)*...
            nanocrystalNumber(sizeIndex)/apertureArea));
    
    end

    % Smear absorbtion spectra with source width
    smearedAbsorbtion = conv(analyticalAbsorbtion, sourceSpectra, 'same');
  
    % Convert back to extinction
    smearedCrossSection = -apertureArea/nanocrystalNumber(sizeIndex)*...
        log(1-smearedAbsorbtion);
    
    figure; 
    subplot(1,2,1); hold on;
    plot(frequencyRange/2/pi, analyticalAbsorbtion, 'b')
    
    plot(frequencyRange/2/pi, smearedAbsorbtion, 'r')
    
    subplot(1,2,2); hold on;
    plot(frequencyRange/2/pi, extinctionCrossSection, 'b')
    
    plot(frequencyRange/2/pi, smearedCrossSection, 'r')
end
