clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    
    
% Actual max mode number is -1
modesToTest = 3;

sizeSteps = 0.5/10^9; % in m

scaleCoreSize = 1; % Fixing core size seems to have quite a large effect on lowest peak

blurAbsorbtion = 0;

% Unkowns are q and Q - firstly determine just using abs. and bandwidth
% Expanded range as padding for convolution
frequencyRange_rad = (25:475)*10^9*2*pi;

%%% Adapted to test on slope range
vLRange = 3200:50:3650;
vTRange = 1500:50:1800;

% Fixed to first three modes
slopeMap = calculateSlopeMap(vLRange, vTRange, modesToTest, 15/10^9, 5/10^9);

% For error fitting
topFraciton = 0.5;

comparisonFreqInds = zeros(length(curveFrequncy),1);

for iFreq = 1:length(curveFrequncy)
    comparisonFreqInds(iFreq) = find(frequencyRange_rad == curveFrequncy(iFreq)*2*pi);
end

% Remove those greater than 300, 2nd peak too complicated for now! 
curveIndsToUse = find(curveFrequncy <= 300*10^9);

testParam = {'systemQ', 'qToUse'};

% For q as quadratic of diameter
% testRange = {2:6:20, 0.8:0.2:1.6};

%%% Add flag to switch fns in future

% For q as fn of volume
testRange = {3:3:15, 0.15:0.05:0.35};

numParamTests = length(testParam);

% Store errors across parameter ranges
errorMapbySize = zeros(numel(slopeMap(:,:,1)), length(testRange{1}), length(testRange{2}),  length(sizesToUse));

extincionCurveBySize = cell(length(sizesToUse), length(testRange{1}), length(testRange{2}));

% Refactored plotting
% Base size x param A - lines param B. For each metric (e.g. extinction) 

colsB = jet(length(testRange{2}));

figure;

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDist = nanocrystalSizeDist{sizeIndex, 1};
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
    [countDist, diameterDist] = interpolatescaleddistribution(countDist, diameterDist, sizeSteps);
    
    coreDist = nanocrystalCore_m(sizeIndex).*diameterDist(:)/ ...
                nanocrystalSize_m(sizeIndex);  
            
    sizeFrequency = countDist/sum(countDist);
    
    % Create blurring function from source
    if freqResolution_Ghz(sizeIndex) == 50
        % convert freq resolution (when FWHM) to sigma
        
        % Seems to work, even though in GHz..?
        sigma = freqResolution_Ghz(sizeIndex)/(2*sqrt(2*log(2)));

        sourceSpectra = 1/(sigma*sqrt(2*pi))*exp(-(((1:200)-100)).^2/...
            (2*(sigma)^2));
        
    else
        error('No conversion for resolution to sigma')
    end
    
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

    % Absorbtion and extinction for bulk    
    [bulkAbsorbtion, bulkExtinctionCrossSection] = calculatesphereabsorbtion(frequencyRange_rad, resonance_rad, ...
        reducedMass, systemQ, qToUse, nanocrystalNumber(sizeIndex), apertureArea);
    
    % Pre determine fixed parameters for each size
    resonanceBySize_rad = zeros(length(diameterDist), modesToTest);
    
    reducedMassBySize = zeros(length(diameterDist), 1);
    
    sphereAreaBySize = 4*pi*(diameterDist/2).^2;
    
    sphereVolumeBySize = 4/3*pi*(diameterDist/2).^3;
      
    for jDiameter = 1:length(diameterDist)

        if scaleCoreSize
            % Scale core diameter given distribution
            
            %%% Assumes core scales as function of shell diameter, 
                % probably not true as each shell will have full range of shell widths
            % Update to full range of core diameters for each shell size
            
            coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDist(jDiameter)/ ...
                nanocrystalSize_m(sizeIndex);
            
            %%% This does actually give same CV for core as shell distribution 
                % not expected...
        else
            coreDiameterScaled_m = nanocrystalCore_m(sizeIndex);
        end
        
        coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

        totalVolume = 4/3*pi*(diameterDist(jDiameter)/2).^3;

        coreMass = coreVolume*CdSeDensity_kgpm3;

        % Take difference in volume for shell mass
        if totalVolume - coreVolume > 0
            shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
        else
           % If core larger than shell, don't let it go negative
           shellMass = 0;
        end
        
        % This will go to zero is shell has zero mass - no absorbtion occurs
        reducedMassBySize(jDiameter) = coreMass*shellMass/(coreMass + shellMass);
        
        % Now using slope map
%         coreFraction = coreMass/(coreMass + shellMass);
% 
%         avgLongVel = CdSeVelocity_mps(1)*coreFraction + ...
%                     CdTeVelocity_mps(1)*(1-coreFraction);
% 
%         avgTransVel = CdSeVelocity_mps(2)*coreFraction + ...
%             CdTeVelocity_mps(2)*(1-coreFraction);
% 
%         % Calculate resonance for this size
%         
%         %%% Set up linear interp across size range based on recipriocal?
%         % This is probably slowest fn.
%         resonanceBySize_rad(jDiameter, :) = calcualtesphereresonance(diameterDist(jDiameter)/2, ...
%                 'sph', 1, modesToTest-1, avgLongVel, avgTransVel, 5*10^9, 10^6, 0)*2*pi;
    end
    
    %%% May be better to make a recursive function for testing more
    %%% parameters, for now, just fix nested  loops
    
    % Do parameter A first
    paramRangeA = testRange{1};

    numRangeTestsA = length(paramRangeA);
    
    measuredExtinctionCurve = ExCrossSecCurve_m2{sizeIndex};
    
    for aRange = 1:numRangeTestsA
        
        % Set up plotting
        subplot(length(sizesToUse), numRangeTestsA, (iSize - 1)*numRangeTestsA + aRange); hold on;

        % Reference cross section
        plot(curveFrequncy/10^9, measuredExtinctionCurve/10^-21, 'm-x')

        normalizationValue = max(measuredExtinctionCurve(curveIndsToUse));
        
        % From bulk
        plot(frequencyRange_rad/2/pi/10^9, bulkExtinctionCrossSection/10^-21, 'k:')

        ylabel('Extinciton')

        xlim([25 475])
        
        if aRange == 1
            title(sprintf('Base size %.1f Testing QF %.1f', nanocrystalSize_m(sizeIndex)*10^9, ...
                paramRangeA(aRange)))
        else
            title(sprintf('Testing QF %.1f', paramRangeA(aRange)))
        end
        
        % Set A param being tested
        switch testParam{1}
            case 'qToUse'
            % Testing charge relationship
            % e as a quadratic fun of diameter in nm
%             qToUse = paramRangeA(aRange)*(diameterDist*10^9).^2* ...
%                 (1.602176634*10^-19);

            % e as linear fun of volume in nm^3
                qToUse = paramRangeA(aRange)*sphereVolumeBySize*10^27*...
                    (1.602176634*10^-19);

            case 'systemQ'
                systemQ = paramRangeA(aRange);

            otherwise

            error('Test param not implemented')
        end    
        
        % Then do parameter B
        paramRangeB = testRange{2};

        numRangeTestsB = length(paramRangeB);

        colsB = jet(numRangeTestsB);
        
        for bRange = 1:numRangeTestsB
            % Set B param being tested
            switch testParam{2}
                case 'qToUse'
                % Testing charge relationship
                % e as a quadratic fun of diameter in nm
%                 qToUse = paramRangeB(bRange)*(diameterDist*10^9).^2* ...
%                     (1.602176634*10^-19);

                    % e as linear fun of volume in nm^3
                    qToUse = paramRangeB(bRange)*sphereVolumeBySize*10^27*...
                        (1.602176634*10^-19);

                case 'systemQ'
                    systemQ = paramRangeB(bRange);

                otherwise
                    error('Test param not implemented') 
            end
            
            % Now calculate across slope distribution
            errorMap = zeros(length(vLRange), length(vTRange));
            
            extinctionCrossSection = zeros(length(frequencyRange_rad), numel(errorMap));
            
            for jSlope = 1:numel(errorMap)
                
                % Was previously broken down by size, but removed for now
                analyticAbsorbtion = zeros(length(frequencyRange_rad),1);
                
                % Do for each mode
                for kMode = 1:modesToTest
                    
                    % Get resonant freq across size dist for this mode
                    tempMap = slopeMap(:, :, kMode);
                    
                    resonanceBySize = 1./(diameterDist/2) * tempMap(jSlope);
                
                    % Could probably vectorize diameter calc rather than using loop
                    
                    % Now calculate for each diameter and sum
                    for lDiameter = 1:length(diameterDist) 
                        
                        % Calculate for each mode and sum
                        if reducedMassBySize(lDiameter) > 0
                            
                            tempAbs = calculatesphereabsorbtion(...
                                frequencyRange_rad, resonanceBySize(lDiameter), reducedMassBySize(lDiameter), systemQ, qToUse(lDiameter), ...
                                nanocrystalNumber(sizeIndex)*sizeFrequency(lDiameter), apertureArea);    

                            analyticAbsorbtion = 1 - (1 - analyticAbsorbtion) .* (1 - tempAbs);
                        end
                    end
                end
                

                % Calc blurred
                if blurAbsorbtion
                    analyticAbsorbtion = conv(analyticAbsorbtion, sourceSpectra, 'same');
                end

                
                % Calculate extinction from absorbtion curve
                extinctionCrossSection(:,jSlope) = -apertureArea/nanocrystalNumber(sizeIndex).*...
                        log(1-analyticAbsorbtion);
                    
                %%% This cheated and lead to points with no comparison    
                % Get top parts of extinction curve
                %topInds = find(extinctionCrossSection(:,jSlope) > max(extinctionCrossSection(:,jSlope))*topFraciton);
                
                % Select those that are points to compare
                %[testInds, measureInds] = intersect(comparisonFreqInds, topInds);
                
                % Store SSE error
                %errorMap(jSlope) = sum((measuredExtinctionCurve(measureInds) - extinctionCrossSection(testInds,jSlope)').^2);
                errorMap(jSlope) = sum((measuredExtinctionCurve(curveIndsToUse) - ... 
                    extinctionCrossSection(comparisonFreqInds(curveIndsToUse),jSlope)').^2);
                
                if imag(errorMap(jSlope)) > 0
                   b = 1; 
                end
            end
            
            % Find lowest error
            [~, minInd] = min(errorMap(:));
            
            [vLInd, vTInd] = ind2sub(size(slopeMap), minInd);
            
            % Plot this and note freqs
            plot(frequencyRange_rad/2/pi/10^9, extinctionCrossSection(:, minInd)/10^-21, 'color', colsB(bRange,:))
            
            [maxEx, maxInd] = max(extinctionCrossSection(:,minInd));
            
%             text(frequencyRange_rad(maxInd)/2/pi/10^9, maxEx/10^-21, ...
%                 sprintf('%i, %i', vLRange(vLInd), vTRange(vTInd)));

            text(frequencyRange_rad(maxInd)/2/pi/10^9, maxEx/10^-21, ...
                sprintf('Charge %.2f', paramRangeB(bRange)));


            % Store for full comp
            errorMapBySize(:, aRange, bRange, iSize) = errorMap(:)/normalizationValue;
            
            extincionCurveBySize{iSize, aRange, bRange} = extinctionCrossSection;
        end
    end
end

combinedErrorMap = sum(errorMapBySize, 4);

[~, minInd] = min(combinedErrorMap(:));
            
[slopeInd, aInd, bInd] = ind2sub(size(combinedErrorMap), minInd);

paramRangeA(aInd)

paramRangeB(bInd)

[vLInd, vTInd] = ind2sub(size(errorMap), slopeInd);

vLRange(vLInd), vTRange(vTInd)

%% Plot best combo
figure;

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);
    
    subplot(1, length(sizesToUse), iSize); hold on
    
    curvesToUse = extincionCurveBySize{iSize, aInd, bInd};
    
    curvesToUse = curvesToUse(:, slopeInd);
    
    measuredExtinctionCurve = ExCrossSecCurve_m2{sizeIndex};
    
    % Reference cross section
    plot(curveFrequncy/10^9, measuredExtinctionCurve/10^-21, 'm-x')

    ylabel('Extinciton')

    xlim([25 475])
    plot(frequencyRange_rad/2/pi/10^9, curvesToUse/10^-21, 'r')
    
    % Find best of this specific size
    
    tempErrorMap = errorMapBySize(:,:,:,iSize);

    [~, tempMinInd] = min(tempErrorMap(:));

    [tempSlopeInd, tempAInd, tempBInd] = ind2sub(size(tempErrorMap), tempMinInd);

    [tempVLInd, tempVTInd] = ind2sub(size(errorMap), tempSlopeInd);

    curvesToUse = extincionCurveBySize{iSize, tempAInd, tempBInd};
    
    curvesToUse = curvesToUse(:, tempSlopeInd);
    
    plot(frequencyRange_rad/2/pi/10^9, curvesToUse/10^-21, 'b')
    
    [maxEx, maxInd] = max(curvesToUse);

    text(frequencyRange_rad(maxInd)/2/pi/10^9, maxEx/10^-21, ...
        sprintf('%i, %f, %i %i', ...
        paramRangeA(tempAInd), paramRangeB(tempBInd), vLRange(tempVLInd), vTRange(tempVTInd)));
end