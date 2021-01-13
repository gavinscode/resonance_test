clear; close all; clc

% Load data
nanosphere_reference

sizesToUse = [1 4 6];

% Unkowns are q and Q - firstly determine just using abs. and bandwidth
frequencyRange = (100:300)*10^9*2*pi;

% Blurring filtesr from source bandwidth

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get calcualted parameter values
        % Note, all frequenciesa are in radians
    resonance = calcualtesphereresonance(nanocrystalSize_m(sizeIndex)/2, ...
            'sph', 0, 4, CdSeVelocity_mps(1), CdSeVelocity_mps(2), 5*10^9, 10^6, 1)*2*pi;
    
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    %%% Q is an unkown, unsure how it should vary with size
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
