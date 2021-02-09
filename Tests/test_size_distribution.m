clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

sizesWithDists = [1 4 6];    

figure; hold on

diameterStep = 0.1;

for iDist = 1:length(sizesWithDists)
    % Plot course distribution
    diameterDist = nanocrystalSizeDist{sizesWithDists(iDist), 1};
    
    countDist = nanocrystalSizeDist{sizesWithDists(iDist), 2};
    
    m = wmean(diameterDist, countDist);
    
    s = std(diameterDist, countDist);
    
    % Weird, histograms match exactly but mean and SD slighly off...
    [m*10^9 s*10^9 s/m]
    
    subplot(1,3,1); hold on
    plot(diameterDist*10^9, countDist, 'o-');
    subplot(1,3,2); hold on
    plot(diameterDist*10^9, countDist/sum(countDist), 'o-');
    
    % Get edges of bins (between centers)
    binWidths = diff(diameterDist);
    
    binEdges = diameterDist(1:end-1) + binWidths/2;
    
    binEdges = [diameterDist(1) - binWidths(1)/2, binEdges, ...
        diameterDist(end) + binWidths(end)/2];
    
    % Do linear interpolation between points 
    interpCounts = interp1(diameterDist*10^9, countDist, ...
       binEdges(1)*10^9:diameterStep:binEdges(end)*10^9, 'linear', 'extrap');

    interpCounts(interpCounts < 0) = 0;
   
    interpDiameters = (binEdges(1)*10^9:diameterStep:binEdges(end)*10^9)/10^9;
    
    subplot(1,3,3); hold on
    plot(interpDiameters*10^9, interpCounts, 'x-');
    
    % Retake count dist using original bins
    interpCountDist = zeros(length(countDist), 1);
    
    for jCount = 1:length(countDist)
       inds = find(interpDiameters > binEdges(jCount) & ...
            interpDiameters <= binEdges(jCount+1));
        
       interpCountDist(jCount) = sum(interpCounts(inds)); 
    end
    
    subplot(1,3,2)
    plot(diameterDist*10^9, interpCountDist/sum(interpCountDist), 'x-');

    % Now rescale counts to fit original bins
    for jCount = 1:length(countDist)
       inds = find(interpDiameters > binEdges(jCount) & ...
            interpDiameters <= binEdges(jCount+1));
        
       ratio = (interpCountDist(jCount)/sum(interpCountDist)) / ...
           (countDist(jCount)/sum(countDist));
        
       interpCounts(inds) = interpCounts(inds)/ratio;
    end
    
    subplot(1,3,3)
    plot(interpDiameters*10^9, interpCounts, '*-');
    
    m = wmean(interpDiameters, interpCounts);
    
    s = std(interpDiameters, interpCounts);
    
    %[m*10^9 s*10^9 s/m]
    
     % Retake count dist using original bins
    interpCountDist = zeros(length(countDist), 1);
    
    for jCount = 1:length(countDist)
       inds = find(interpDiameters > binEdges(jCount) & ...
            interpDiameters <= binEdges(jCount+1));
        
       interpCountDist(jCount) = sum(interpCounts(inds)); 
    end
    
    subplot(1,3,2)
    plot(diameterDist*10^9, interpCountDist/sum(interpCountDist), '*-');
    
    % Test interpolation fn
%     [newCounts, newBins] = interpolatescaleddistribution(countDist,...
%         diameterDist, diameterStep/10^9);
%     
%     subplot(1,3,3)
%     plot(newBins*10^9, newCounts, 'd-');
end

%% Look at core ratio distributions
% Very funky, not linear

figure; hold on

for iDist = 1:length(sizesWithDists)
    sizeIndex = sizesWithDists(iDist);
    
    % Plot course distribution
    diameterDist = nanocrystalSizeDist{sizeIndex, 1};
    
    countDist = nanocrystalSizeDist{sizeIndex, 2};
    
    coreDiameterScaled_m = nanocrystalCore_m(sizeIndex).*diameterDist/ ...
            nanocrystalSize_m(sizeIndex);
        
    coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

    totalVolume = 4/3*pi*(diameterDist/2).^3;
    
    volumeRatio = zeros(length(diameterDist), length(diameterDist));
    
    volumeRatioCounts = zeros(length(diameterDist), length(diameterDist));
    
    for jRatio = 1:length(diameterDist)
        volumeRatio(jRatio,:) = coreVolume(jRatio)./totalVolume;
        
        volumeRatioCounts(jRatio,:) = countDist(jRatio).*countDist/sum(countDist);
    end
    
    [volumeRatio, inds] = sort(volumeRatio(:));
    
    volumeRatio(volumeRatio > 1) = 1;
    
    plot(volumeRatio, volumeRatioCounts(inds))
    
    %plot(coreVolume./totalVolume, countDist);
end