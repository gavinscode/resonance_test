clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    

interpolateSizeDist = 0; % If not set, uses original distrabution

sizeSteps = 1/10^9; % in m

% typeOfSolution = '3Coeffs';

typeOfSolution = 'IntensityPerSize';

frequencyRange_rad = (50:400)*10^9*2*pi;

figure;

for iSize = 3; 1:length(sizesToUse);
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDistOrig = nanocrystalSizeDist{sizeIndex, 1};
    
    countDistOrig = nanocrystalSizeDist{sizeIndex, 2};
    
    if interpolateSizeDist
        [countDist, diameterDist] = interpolatescaleddistribution(countDistOrig, diameterDistOrig, sizeSteps);
    else
        diameterDist = diameterDistOrig;

        countDist = countDistOrig;
    end
    
    sizeFrequencyOrig = countDistOrig/sum(countDistOrig);

    sizeFrequency = countDist/sum(countDist);
    
    coreDist = nanocrystalCore_m(sizeIndex).*diameterDist(:)/ ...
                nanocrystalSize_m(sizeIndex);  
    
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
        reducedMass, systemQ, qToUse, nanocrystalNumber(sizeIndex), apertureArea, 1, []);
    
    testPower = calculateresonatorpower(frequencyRange_rad, resonance_rad, qToUse/sqrt(reducedMass), ...
        systemQ);
        
    % Assumes drive power is 1
    testExtinction = testPower/(0.5*VACCUM_PERMITIVITY*LIGHT_SPEED);
    
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
    
    % Calculate extinction for sphere mixture
    [~, simulatedExtinctionCrossSection] = calculatespheremixtureabsorbtion(frequencyRange_rad, resonanceBySize_rad, ...
        reducedMassBySize, systemQ, qToUse, nanocrystalNumber(sizeIndex)*sizeFrequency, apertureArea, 1, []);
    
    figure;
    plot(curveFrequncy/10^9, measuredExtinctionCurve/10^-21, 'm-x'); hold on
    
    plot(frequencyRange_rad/2/pi/10^9, bulkExtinctionCrossSection/10^-21, 'g')

    plot(frequencyRange_rad/2/pi/10^9, testExtinction/10^-21, 'r--')
    
    plot(frequencyRange_rad/2/pi/10^9, simulatedExtinctionCrossSection/10^-21, 'g')
    
    switch typeOfSolution
    
        % Use solver for 3 coefficients
        case '3Coeffs'
            % Parameters rescalled to similar ranges
            x0 = [resonantSlope/1000, qToUse/sqrt(reducedMass)*10^6, systemQ];

            f = @(x)(calculateerror_equalcoefficients(x, simulatedExtinctionCrossSection, diameterDistOrig, ...
                nanocrystalNumber(sizeIndex)*sizeFrequencyOrig, apertureArea, frequencyRange_rad, 10^21));
    
        % Use solver allowing for varying intensity    
        case 'IntensityPerSize'
        %%% Soln values quite sensitive to initial parameter choices
            % Acceptable space deffinetely worth considering...
        x0 = [resonantSlope/1000, qToUse/sqrt(reducedMass)*10^6*ones(1,length(diameterDistOrig)),...
            systemQ];

        f = @(x)(calculateerror_varyingintensity(x, simulatedExtinctionCrossSection, diameterDistOrig, ...
            nanocrystalNumber(sizeIndex)*sizeFrequencyOrig, apertureArea, frequencyRange_rad, 10^21));
    end
    
    options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10);
    
    % Use lower bound at zero for all
    solution = lsqnonlin(f, x0, zeros(length(x0),1), [], options)
    
    solution - x0
    
    solution = [solution(1)*1000 solution(2:end-1)/10^6 solution(end)];
    
    %%% Make a function that does this (take from solved function as well)
    resonance_rad = 1./(diameterDistOrig/2) * solution(1);
    
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
        solution(2:end-1), solution(end), nanocrystalNumber(sizeIndex)*sizeFrequencyOrig, apertureArea, 1, []);
    
    plot(frequencyRange_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'r--')
end

% Creates extinction curve from resonator. Multiple sizes, but all use equal coefficients
function sse = calculateerror_equalcoefficients(x, extinctionRef, diameters, number, area, ...
    frequencyRange, conversion)
  
    % For solver
    % Expects, x(1): slope, x(2): intensity (q/mass), x(3): gamma
    % Expects all frequency values in radians and diameters in m
    
    if length(x) ~= 3
       error('wrong number of test parameters') 
    end
    
    x = [x(1)*1000 x(2)/10^6 x(3)];
    
    resonance_rad = 1./(diameters/2) * x(1);
    
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRange, resonance_rad, ...
        x(2), x(3), number, area, 1, []);

    
    % Apply conversion to error to keep it above solver tolerance.
    
    % Top for fmincon - requires single value
%     sse = sum((conversion*(extinctionCrossSection-extinctionRef)).^2);

    % Bottom for lsqnonlin - requires array of differences, then matches fmincon
    sse = conversion*(extinctionCrossSection-extinctionRef);
end

% Creates extinction curve from resonator. Multiple sizes and varying intensity - x(2)
function sse = calculateerror_varyingintensity(x, extinctionRef, diameters, number, area, ...
    frequencyRange, conversion)

    % For solver
    % Expects, x(1): slope, x(2): intensity (q/mass), x(3): gamma
    % Expects all frequency values in radians and diameters in m
    
    if length(x) ~= length(diameters)+2
       error('wrong number of test parameters') 
    end
    
    x = [x(1)*1000 x(2:end-1)/10^6 x(end)];
    
    resonance_rad = 1./(diameters/2) * x(1);
    
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRange, resonance_rad, ...
        x(2:end-1), x(end), number, area, 1, []);
    
    % Apply conversion to error to keep it above solver tolerance.
    
    % Top for fmincon - requires single value
%     sse = sum((conversion*(extinctionCrossSection-extinctionRef)).^2);

    % Bottom for lsqnonlin - requires array of differences, then matches fmincon
    sse = conversion*(extinctionCrossSection-extinctionRef);
end

