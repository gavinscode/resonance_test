clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    
    
sizeSteps = 1/10^9; % in m

% Unkowns are q and Q - firstly determine just using abs. and bandwidth
% Expanded range as padding for convolution
frequencyRange_rad = (25:475)*10^9*2*pi;


figure;

for iSize = 3; 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDist = nanocrystalSizeDist{sizeIndex, 1};
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
%     [countDist, diameterDist] = interpolatescaleddistribution(countDist, diameterDist, sizeSteps);
%     
%     coreDist = nanocrystalCore_m(sizeIndex).*diameterDist(:)/ ...
%                 nanocrystalSize_m(sizeIndex);  
%             
    sizeFrequency = countDist/sum(countDist);
    
    measuredExtinctionCurve = ExCrossSecCurve_m2{sizeIndex};
    
    % Calculate parameters from bulk
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    %%% For now, use Q and q from bulk
    resonance_rad = nanocrystalFreqResonance_hz(sizeIndex)*2*pi;
    
    % Q is an unknown, unsure how it should vary with size
    systemQ = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);

    systemSpring = resonance_rad^2*reducedMass;

    systemDamp = resonance_rad*reducedMass/systemQ;

    % q is an unknown and should vary with the square of size 
    qToUse = sqrt(nanocrystalThetaEx_m2(sizeIndex) * resonance_rad * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED/systemQ);  

    % Get slope of main size
    resonantSlope = resonance_rad/(1/(nanocrystalSize_m(sizeIndex)/2));
    
    % Absorbtion and extinction for bulk    
    [bulkAbsorbtion, bulkExtinctionCrossSection, bulkPower] = calculatesphereabsorbtion(frequencyRange_rad, resonance_rad, ...
        reducedMass, systemQ, qToUse, nanocrystalNumber(sizeIndex), apertureArea, 1);
    
    % Make coarse synthetic and plot against original and single
    %%% Fixed Q and q - should be easy to fit...
    
    resonanceBySize_rad = zeros(length(diameterDist), 1);
    
    reducedMassBySize = zeros(length(diameterDist), 1);
    
    for jDiameter = 1:length(diameterDist)

        coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDist(jDiameter)/ ...
            nanocrystalSize_m(sizeIndex);
        
        coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

        totalVolume = 4/3*pi*(diameterDist(jDiameter)/2).^3;

        coreMass = coreVolume*CdSeDensity_kgpm3;

        shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
        
        % This will go to zero is shell has zero mass - no absorbtion occurs
        reducedMassBySize(jDiameter) = coreMass*shellMass/(coreMass + shellMass);
        
        % Get resonant frequncy
        resonanceBySize_rad(jDiameter) = 1./(diameterDist(jDiameter)/2) * resonantSlope;
    end
    
    analyticAbsorbtion = zeros(length(frequencyRange_rad),1);
    
    for jDiameter = 1:length(diameterDist) 
                       
        tempAbs = calculatesphereabsorbtion(...
            frequencyRange_rad, resonanceBySize_rad(jDiameter), reducedMassBySize(jDiameter), systemQ, qToUse, ...
            nanocrystalNumber(sizeIndex)*sizeFrequency(jDiameter), apertureArea, 1);    

        analyticAbsorbtion = 1 - (1 - analyticAbsorbtion) .* (1 - tempAbs);
    end
      
    extinctionCrossSection = -apertureArea/nanocrystalNumber(sizeIndex).*...
        log(1-analyticAbsorbtion);
    
    figure;
    plot(curveFrequncy/10^9, measuredExtinctionCurve/10^-21, 'm-x'); hold on
    
    plot(frequencyRange_rad/2/pi/10^9, bulkExtinctionCrossSection/10^-21, 'k:')

    plot(frequencyRange_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'g')
    
    testPower = calculateresonatorpower(frequencyRange_rad, resonance_rad, 1*qToUse/sqrt(reducedMass), nanocrystalFreqBandwidth_hz(sizeIndex)*2*pi);

    % testExtinction
    
    %To sum: testAbsorbtion
    
    %To fit: testExtinction
    
    % Solve all at once or in steps
    
end

