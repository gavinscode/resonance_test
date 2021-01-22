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