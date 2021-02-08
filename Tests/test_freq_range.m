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

%% Make map of curve slopes
%%% Note slopes are from reciprical of radius
% Note, this is for bulk sphere, no core mass averaging

testLong = 3200:50:4000;

testTrans = 1200:50:2000;

% Get resonant freq slop range of speed combinations
slopeMap = zeros(length(testLong), length(testTrans), modesToTest);

for jLong = 1:length(testLong)

    for kTrans = 1:length(testTrans)

        % Test two sizes to get slope (sizes aren't really important)
        highFreqs = calcualtesphereresonance(nanocrystalSize_m(sizesToUse(1))/2, ...
            'sph', 1, modesToTest-1, testLong(jLong), testTrans(kTrans), 5*10^9, 10^6, 0)*2*pi;
        
        lowFreqs = calcualtesphereresonance(nanocrystalSize_m(sizesToUse(end))/2, ...
            'sph', 1, modesToTest-1, testLong(jLong), testTrans(kTrans), 5*10^9, 10^6, 0)*2*pi;
        
        % Calc slope
        xTravel = 1/(nanocrystalSize_m(sizesToUse(1))/2) - 1/(nanocrystalSize_m(sizesToUse(end))/2);
        
        yTavel = highFreqs - lowFreqs;
        
        slopeMap(jLong, kTrans, :) = yTavel/xTravel;
    end
end  

% slopes = calculateSlopeMap(testLong, testTrans, modesToTest, nanocrystalSize_m(sizesToUse(end)), nanocrystalSize_m(sizesToUse(1)));
% 
% sum(slopeMap(:) - slopes(:))

mode1 = slopeMap(:,:,1);
mode2 = slopeMap(:,:,2);
mode3 = slopeMap(:,:,3);

figure; hold on;
[n ,x] = hist(mode2(:)./mode1(:));
plot(x , n)
[n ,x] = hist(mode3(:)./mode1(:));
plot(x , n)

[n ,x] = hist(mode3(:)./mode2(:));
plot(x , n, ':')

nanocrystal2ndFreqResonance_hz./nanocrystalFreqResonance_hz

nanocrystal3rdFreqResonance_hz./nanocrystalFreqResonance_hz

nanocrystal3rdFreqResonance_hz./nanocrystal2ndFreqResonance_hz

%%
warning('Adjusted freqs on 1st and 4th')
nanocrystalFreqResonance_hz = [230 241 175 200 170 165]*10^9;

figure;

sizeCols = jet(6);

for iMode = 1:modesToTest
    subplot(2,modesToTest,iMode)
    
    imshow(slopeMap(:,:,iMode)/max(slopeMap(:)))
    hold on
    
    axis on
    set(gca, 'YTickLabel', testTrans, 'YTick', 1:length(testTrans), ...
        'XTickLabel', testLong, 'XTick', 1:length(testLong))
    
    subplot(2,modesToTest,iMode+modesToTest); hold on
    
    testSizes = (5:15)/10^9;
    
    for jLong = 1:length(testLong)

        for kTrans = 1:length(testTrans)
            freqs = 1./(testSizes/2) * slopeMap(jLong, kTrans, iMode);
            
            plot(1./(testSizes*10^9), freqs/2/pi/10^9)
        end
    end
    
    xlim([1/15 1/6])
    title('reciprical diameter')
    ylim([100 400])
    
    switch iMode
        case 1
            freqsToUse = nanocrystalFreqResonance_hz;
        
        case 2
            freqsToUse = nanocrystal2ndFreqResonance_hz;    
            
        case 3
            freqsToUse = nanocrystal3rdFreqResonance_hz;
    end
    
    plot(1./(nanocrystalSize_m*10^9), freqsToUse/10^9, 'xk')
    plot(1./(nanocrystalSize_m(sizesToUse)*10^9), freqsToUse(sizesToUse)/10^9, 'ok')

    subplot(2,modesToTest,iMode);   
        
    for jSize = 1:length(nanocrystalSize_m)
        if ~isnan(freqsToUse(jSize))
            slope = (freqsToUse(jSize)*2*pi)./(1/(nanocrystalSize_m(jSize)/2));

            temp = slopeMap(:, :, iMode);

            closeInds = find(abs(temp(:)-slope) < 100);
            
            [xInd, yInd] = ind2sub(size(slopeMap(:,:,iMode)), closeInds);

            plot(yInd+0.1*jSize-0.3, xInd+0.1*jSize-0.3, '.', 'color', sizeCols(jSize,:));
        end
    end
    
end