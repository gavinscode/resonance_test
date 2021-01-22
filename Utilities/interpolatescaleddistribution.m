function [countsOut, binsOut] = interpolatescaleddistribution(countsIn, binsIn, newStep)
    % Bins in and out represent centers

    % Get edges of bins (between centers)
    binWidths = diff(binsIn);
    
    binEdges = binsIn(1:end-1) + binWidths/2;
    
    binEdges = [binsIn(1) - binWidths(1)/2, binEdges, ...
        binsIn(end) + binWidths(end)/2];
    
    % Do linear interpolation between end points of bin edges
    countsOut = interp1(binsIn, countsIn, ...
       min(binEdges):newStep:max(binEdges), 'linear', 'extrap');

    binsOut = min(binEdges):newStep:max(binEdges);

    % Catch in case interpolation goes below zero
    binsOut(countsOut < 0) = [];
    
    countsOut(countsOut < 0) = [];
    
    % Get dist from inpolant using original bins
    interpCountDist = zeros(length(countsIn), 1);
    
    for jCount = 1:length(countsIn)
       inds = find(binsOut > binEdges(jCount) & ...
            binsOut <= binEdges(jCount+1));
        
       interpCountDist(jCount) = sum(countsOut(inds)); 
    end
    
    % Now rescale counts to fit original bins
    for jCount = 1:length(countsIn)
       inds = find(binsOut > binEdges(jCount) & ...
            binsOut <= binEdges(jCount+1));
        
       ratio = (interpCountDist(jCount)/sum(interpCountDist)) / ...
           (countsIn(jCount)/sum(countsIn));
        
       countsOut(inds) = countsOut(inds)/ratio;
    end
end

