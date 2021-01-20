clear; close all; clc

% Load data
nanosphere_reference

absorbtion_reference

sizesWithDists = [1 4 6];    

figure; hold on

for iDist = 1:length(sizesWithDists)
    diameterDist = nanocrystalSizeDist{sizesWithDists(iDist), 1};
    
    countDist = nanocrystalSizeDist{sizesWithDists(iDist), 2};
    
    plot(diameterDist*10^9, countDist/sum(countDist));
end

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = 1; [1 4 6];    
    
% Unkowns are q and Q - firstly determine just using abs. and bandwidth
% Expanded range as padding for convolution
frequencyRange_rad = (50:450)*10^9*2*pi;

% testParam = {'qToUse'};
% 
% testRange = {1:0.2:2};

testParam = {'systemQ', 'qToUse'};

testRange = {2:2:10, 0.8:0.2:1.6};

testRange = {[2 20], [0.8 1.6]};

numParamTests = length(testParam);
    
for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDist = nanocrystalSizeDist{sizeIndex, 1}*10^9;
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
    % Add interpolation function
    countDist = interp1(diameterDist, countDist, ...
        min(diameterDist):0.5:max(diameterDist), 'linear', 0);
    
    %%% Check this agrees with histogram predictions...
    diameterDist = (min(diameterDist):0.5:max(diameterDist))/10^9;
    
    sizeFrequency = countDist/sum(countDist);
    
    % Create blurring function from source
    if freqResolution_Ghz(sizeIndex) == 50
        % convert freq resolution (when FWHM) to sigma
        sigma = freqResolution_Ghz(sizeIndex)/(2*sqrt(2*log(2)));

        sourceSpectra = 1/(sigma*sqrt(2*pi))*exp(-((1:200)-100).^2/...
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
    resonanceBySize_rad = zeros(length(diameterDist), 1);
    
    reducedMassBySize = zeros(length(diameterDist), 1);
    
    for jDiameter = 1:length(diameterDist)

        % Scale core diameter given distribution
        %%% Note, assumse core directly scales with total diamter, may not be true
        coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDist(jDiameter)/ ...
            nanocrystalSize_m(sizeIndex);

        coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

        totalVolume = 4/3*pi*(diameterDist(jDiameter)/2).^3;

        coreMass = coreVolume*CdSeDensity_kgpm3;

        % Take difference in volume
        shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;

        reducedMassBySize(jDiameter) = coreMass*shellMass/(coreMass + shellMass);

        coreFraction = coreMass/(coreMass + shellMass);

        avgLongVel = CdSeVelocity_mps(1)*coreFraction + ...
                    CdTeVelocity_mps(1)*(1-coreFraction);

        avgTransVel = CdSeVelocity_mps(2)*coreFraction + ...
            CdTeVelocity_mps(2)*(1-coreFraction);
        
        % Calculate resonance for this size
        resonanceBySize_rad(jDiameter) = calcualtesphereresonance(diameterDist(jDiameter)/2, ...
                'sph', 1, 0, avgLongVel, avgTransVel, 5*10^9, 10^6, 0)*2*pi;
    end
    
    %%% May be better to make a recursive function for testing more
    %%% parameters, for now, just fix nested loops
    
    % Do parameter A first
    paramRangeA = testRange{1};

    numRangeTestsA = length(paramRangeA);
        
    for aRange = 1:numRangeTestsA

        figure;
        
        % Then do parameter B
        paramRangeB = testRange{2};

        numRangeTestsB = length(paramRangeB);

        for bRange = 1:numRangeTestsB
            % Now calculate across diameter distribution
            analyticAbsorbtion = zeros(length(diameterDist), length(frequencyRange_rad));

            analyticExtinctionCrossSection = zeros(length(diameterDist), length(frequencyRange_rad));

            for jDiameter = 1:length(diameterDist)
                
                % Set A param being tested
                switch testParam{1}
                    case 'qToUse'
                    % Testing charge relationship
                    % e as a quadratic fun of diameter in nm
                    qToUse = paramRangeA(aRange)*(diameterDist(jDiameter)*10^9)^2* ...
                        (1.602176634*10^-19);

                    case 'systemQ'
                        systemQ = paramRangeA(aRange);

                    otherwise
                        
                    error('Test param not implemented')
                end
                
                % Set B param being tested
                switch testParam{2}
                    case 'qToUse'
                    % Testing charge relationship
                    % e as a quadratic fun of diameter in nm
                    qToUse = paramRangeB(bRange)*(diameterDist(jDiameter)*10^9)^2* ...
                        (1.602176634*10^-19);

                    case 'systemQ'
                        systemQ = paramRangeB(bRange);

                    otherwise
                        error('Test param not implemented') 
                end
                
                [analyticAbsorbtion(jDiameter, :), analyticExtinctionCrossSection(jDiameter, :)] = calculatesphereabsorbtion(...
                    frequencyRange_rad, resonanceBySize_rad(jDiameter), reducedMassBySize(jDiameter), systemQ, qToUse, ...
                    nanocrystalNumber(sizeIndex)*sizeFrequency(jDiameter), apertureArea);    

                subplot(numRangeTestsB, 2, bRange*2-1); hold on;
                plot(frequencyRange_rad/2/pi/10^9, analyticAbsorbtion(jDiameter, :)*100)

                %subplot(numRangeTestsB, 2, bRange*2); hold on;
                %plot(frequencyRange_rad/2/pi/10^9, analyticExtinctionCrossSection(jDiameter, :)/10^-21)
            end

            % Plot combined absorbtion
            subplot(numRangeTestsB,2,bRange*2-1); hold on;
            totalAbsorbtion = sum(analyticAbsorbtion);

            plot(frequencyRange_rad/2/pi/10^9, totalAbsorbtion*100, 'r')

            % From bulk
            plot(frequencyRange_rad/2/pi/10^9, bulkAbsorbtion*100, 'g')

            % Add measured absorbtion curve
            if sizeIndex == 6
                % Use proved intensity for 13 nm
                directAbs = 1 - I_trans_13_mv./I_source_mv;

                plot(curveFrequncy/10^9, directAbs*100, 'b-x')
            else
               % Otherise calculate intensity 
               intensityCalc = exp(ExCrossSecCurve_m2{sizeIndex}*-nanocrystalNumber(sizeIndex)/apertureArea) .* ...
                   I_source_mv;

               absCalc = 1 - intensityCalc./I_source_mv;

               plot(curveFrequncy/10^9, absCalc*100, 'b-x')
            end

            % Calc and plot blurred absorbtion
            blurredAbsorbtion = conv(totalAbsorbtion, sourceSpectra, 'same');

            plot(frequencyRange_rad/2/pi/10^9, blurredAbsorbtion*100, 'm')

            title(sprintf('Testing %s %.1f', testParam{2}, paramRangeB(bRange)))

            ylabel('Abosorbtion')

            xlim([100 400])

            % Now plot extinction cross section
            subplot(numRangeTestsB,2,bRange*2); hold on;

            % Reference cross section
            plot(curveFrequncy/10^9, ExCrossSecCurve_m2{sizeIndex}/10^-21, 'b-x')

            % From bulk
            plot(frequencyRange_rad/2/pi/10^9, bulkExtinctionCrossSection/10^-21, 'g')

            % Cross section from combined absorbtion curve
            totalExtinctionCrossSection = -apertureArea/nanocrystalNumber(sizeIndex).*...
                log(1-totalAbsorbtion);

            plot(frequencyRange_rad/2/pi/10^9, totalExtinctionCrossSection/10^-21, 'r')

            % Calculate and plot extinction from blurred absorbtion curve;
            blurredExtinction = -apertureArea/nanocrystalNumber(sizeIndex).*...
                log(1-blurredAbsorbtion);

            plot(frequencyRange_rad/2/pi/10^9, blurredExtinction/10^-21, 'm')

            ylabel('Extinciton')

            xlim([100 400])
        end
    
        subplot(numRangeTestsB,2,2);

        title(sprintf('Base size %.1f Testing %s %.1f', nanocrystalSize_m(sizeIndex)*10^9, ...
            testParam{1}, paramRangeA(aRange)))
    end
end
