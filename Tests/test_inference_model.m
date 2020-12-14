clear; close all; clc

% Load data
nanosphere_reference

sizesToUse = [1 4 6];

% Unkowns are q and Q - firstly determine just using abs. and bandwidth
frequencyRange = (100:300)*10^9*2*pi;

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get calcualted parameter values
        % Note, all frequenciesa are in radians
    resonance = calcualtesphereresonance(nanocrystalSize_m(sizeIndex)/2, ...
            'sph', 1, 2, CdSeVelocity_mps(1), CdSeVelocity_mps(2), 100*10^9)*2*pi;
    
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    systemQ = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);
    
    systemSpring = resonance^2*reducedMass;
    
    systemDamp = resonance*reducedMass/systemQ;

    %%% Add solution for Q    
    qToUse = sqrt(nanocrystalThetaEx_m2(sizeIndex) * nanocrystalFreqResonance_hz(sizeIndex) * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED/systemQ);  
    
    analyticalAbsorbtion = zeros(length(testFrequncies_hz), 1);
    
    for kFreq = 1:length(testFrequncies_hz)
        % Eqn 7
        analyticalAmplitude = qToUse./(reducedMass*...
            sqrt((resonance^2 - frequencyRange(jFreq)^2)^2 + ...
            (resonance(jDiam)*frequencyRange(jFreq)/systemQ).^2));

        % Eqn 9
        analyticalPower = resonance * ...
            frequencyRange(jFreq)^2 * reducedMass * ...
            analyticalAmplitude^2 / (2*systemQ);

        powerFlux = 0.5*VACCUM_PERMITIVITY*LIGHT_SPEED;

        % Eqn 10
        absorbtionCrossSection = analyticalPower/powerFlux;
        
        %%% Convert to absorbtion - sum absorbtion across pulse
            % Convert back to extinciton
        
       % Eqn 13 - solution for absorbtion
       analyticalAbsorbtion(jFreq) = (1-exp(-absorbtionCrossSection*...
            nanocrystalNumber/excitationArea));
    
    end

   %%% Convoluve with 50 GHz bandwidth here...
   
end