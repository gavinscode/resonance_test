clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    

% Actual max mode number is -1
modesToTest = 3;

modeCols = lines(modesToTest);

modeScale = [0.8 1 1.2];

sizeSteps = 0.5/10^9; % in m

scaleCoreSize = 1; % Fixing core size seems to have quite a large effect on lowest peak

frequencyRange_rad = (50:450)*10^9*2*pi;

% Plot ranges on existing curves

figure;

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDist = nanocrystalSizeDist{sizeIndex, 1};
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
    % Don't interpoalte for now
    %[countDist, diameterDist] = interpolatescaleddistribution(countDist, diameterDist, sizeSteps);
    
    % Just use smallest, most common, and largest in dist
    [~, minD] = min(diameterDist);
    
    [~, maxD] = max(diameterDist);
    
    [~, modeD] = max(countDist);
    
    diameterDist = diameterDist([minD, modeD, maxD]);
    
    countDist = countDist([minD, modeD, maxD]);
    
    areaDist = 4*pi*(diameterDist/2).^2;
    
    volumeDist = 4/3*pi*(diameterDist/2).^3;
    
    % Pre determine fixed parameters for each size
    resonanceBySize_rad = zeros(length(diameterDist), modesToTest);
    
    for jDiameter = 1:length(diameterDist)

        if scaleCoreSize
            % Scale core diameter given distribution
            %%% Note, assumse core directly scales with total diamter, may not be true
            coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDist(jDiameter)/ ...
                nanocrystalSize_m(sizeIndex);
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
        
        coreFraction = coreMass/(coreMass + shellMass);

        avgLongVel = CdSeVelocity_mps(1)*coreFraction + ...
                    CdTeVelocity_mps(1)*(1-coreFraction);

        avgTransVel = CdSeVelocity_mps(2)*coreFraction + ...
            CdTeVelocity_mps(2)*(1-coreFraction);
        
        % Calculate resonance for this size
        
        %%% Set up linear interp across size range based on recipriocal?
        % This is probably slowest fn.
        resonanceBySize_rad(jDiameter, :) = calcualtesphereresonance(diameterDist(jDiameter)/2, ...
                'sph', 1, modesToTest-1, avgLongVel, avgTransVel, 5*10^9, 10^6, 0)*2*pi;
    end
    
    % Plot curve and ranges for each mode
    
    subplot(1, length(sizesToUse), iSize); hold on
    
    plot(curveFrequncy/10^9, ExCrossSecCurve_m2{sizeIndex}/10^-21, 'm-x')
    
    meanCurve = mean(ExCrossSecCurve_m2{sizeIndex}/10^-21);
    
    ylabel('Extinciton')

    title(sprintf('Base size %.1f', nanocrystalSize_m(sizeIndex)*10^9));
    
    xlim([50 450])

    for jMode = 1:modesToTest

        plot(resonanceBySize_rad(:,jMode)/2/pi/10^9, ones(3,1)*meanCurve*modeScale(jMode), ...
            'o-', 'color', modeCols(jMode,:));
        
        % Test area weighted
        weightDist = countDist.*areaDist;
        
        % Take weighted mean between each extreme
        highMean = wmean(resonanceBySize_rad(2:3,jMode), weightDist(2:3)');
        
        plot(highMean/2/pi/10^9, ones(2,1)*meanCurve*modeScale(jMode), ...
            'x', 'color', modeCols(jMode,:));
        
        lowMean = wmean(resonanceBySize_rad(1:2,jMode), weightDist(1:2)');
        
        plot(lowMean/2/pi/10^9, ones(2,1)*meanCurve*modeScale(jMode), ...
            'x', 'color', modeCols(jMode,:));
        
        % Test volume weighted
        weightDist = countDist.*volumeDist;
        
        % Take weighted mean between each extreme
        highMean = wmean(resonanceBySize_rad(2:3,jMode), weightDist(2:3)');
        
        plot(highMean/2/pi/10^9, ones(2,1)*meanCurve*modeScale(jMode), ...
            'd', 'color', modeCols(jMode,:));
        
        lowMean = wmean(resonanceBySize_rad(1:2,jMode), weightDist(1:2)');
        
        plot(lowMean/2/pi/10^9, ones(2,1)*meanCurve*modeScale(jMode), ...
            'd', 'color', modeCols(jMode,:));
        
    end
end

% Ideal low freq range from maps
% For 8, 10.4 and 13 nm spheres (1 4 6)
peakRange = [150 250; 150 250; 125 225];

peakRange(:,:,2) = [350 400; 200 350; 175 275];

peakRange(:,:,3) = [NaN NaN; 300 450; 325 400];

% Make map of test values
% Note, this is for bulk sphere, no core mass averaging

testLong = 3000:200:4000;

testTrans = 1000:200:2000;

figure;

overLapMap = zeros(length(testLong), length(testTrans));

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDist = nanocrystalSizeDist{sizeIndex, 1};
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
    % Don't interpolate for now
    %[countDist, diameterDist] = interpolatescaleddistribution(countDist, diameterDist, sizeSteps);
    
    % Just use smallest, most common, and largest in dist
    [~, minD] = min(diameterDist);
    
    [~, maxD] = max(diameterDist);
    
    [~, modeD] = max(countDist);
    
    %Leave out min d, it's never in the range
    diameterDist = mean(diameterDist([modeD maxD]));
    
%     countDist = countDist([modeD]);
%     
%     areaDist = 4*pi*(diameterDist/2).^2;
%     
%     volumeDist = 4/3*pi*(diameterDist/2).^3;
%     
%     weightDist = countDist.*volumeDist;
%     
%     weightDist = weightDist / sum(weightDist);
    
    %%% Just have test first mode
    
    %for aDiameter = 1:length(diameterDist)
    
        % Get resonant freq across range of speed combinations for each size
        resultsMap = zeros(length(testLong), length(testTrans), modesToTest);

        for jLong = 1:length(testLong)

            for kTrans = 1:length(testTrans)

                resultsMap(jLong, kTrans,:) = calcualtesphereresonance(diameterDist/2, ...
                    'sph', 1, modesToTest-1, testLong(jLong), testTrans(kTrans), 5*10^9, 10^6, 0)*2*pi;
            end
        end  

        for jMode = 1:modesToTest
            subplot(modesToTest, length(sizesToUse), (jMode-1)*length(sizesToUse)+iSize); 

            resultsMap_ghz = resultsMap(:,:,jMode)/2/pi/10^9;

            resultsMapScaled = (resultsMap_ghz-100)/400;

            imshow(resultsMapScaled)

            hold on

            title(sprintf('Base size %.1f', nanocrystalSize_m(sizeIndex)*10^9));

            % label points in range
            inds = find(resultsMap_ghz > peakRange(iSize,1,jMode) & resultsMap_ghz < peakRange(iSize,2,jMode));

            %overLapMap(inds) = overLapMap(inds) + weightDist(aDiameter);

            [pointX, pointY] = ind2sub(size(resultsMap), inds);

            plot(pointY, pointX, 'rx');

            axis on
            set(gca, 'XTickLabel', testLong(1:2:end), 'XTick', 1:2:length(testLong))
            set(gca, 'YTickLabel', testTrans(1:2:end), 'YTick', 1:2:length(testTrans))
        end

    %end
end

colorH = colorbar;

set(colorH, 'XTick', 0:0.25:1, 'XTickLabel', 100:100:400)

% Plot overlapping points.
% inds = find(overLapMap == 3);
% 
% [pointX, pointY] = ind2sub(size(resultsMap), inds);
% 
% for iSize = 1:length(sizesToUse)
%     subplot(1, length(sizesToUse), iSize); 
% 
%     plot(pointY, pointX, 'go');
%     
%     
% end

%figure; imshow(overLapMap/6)