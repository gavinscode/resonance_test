function [sse, freqSizeStrucutre, result] = inactivationError(weightsVector, freqSizeStrucutre, sizeWeighting, ...
    inactivationMap, weightsIndexes, powerReference)

    freqSizeStrucutre(weightsIndexes) = weightsVector;

    result = zeros(size(inactivationMap));
    
    for i = 1:size(result,1)
        
        for j = 1:size(result,2)
            inds = find(freqSizeStrucutre(i,:) <= powerReference(j) & freqSizeStrucutre(i,:) > 0);
            
            result(i,j) = sum(sizeWeighting(inds));
        end
    end
    
    sse = sum((result(:)-inactivationMap(:)).^2); %sum
end

