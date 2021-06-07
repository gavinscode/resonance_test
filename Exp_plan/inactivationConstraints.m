function [c, ceq] = inactivationConstraints(weightsVector, freqSizeStrucutre, weightsIndexes,...
        minimaType)

    freqSizeStrucutre(weightsIndexes) = weightsVector;

    % change to linear power for convience
    if any(isinf(freqSizeStrucutre(:)))
        freqSizeStrucutre = 10.^(freqSizeStrucutre);
    end
    
    freqSizeStrucutre = ceil(freqSizeStrucutre);
    
    c = [];
    
    % c is inequality constraint, solves for c(x) < 0
    % calculate c for each size
    
    %%% Set log10/integer on equality constraint
    % probably works best with integer points...
    
%     uniThresholds = unique(freqSizeStrucutre(freqSizeStrucutre(:) > 0));
%     
%     c = zeros(size(freqSizeStrucutre,2),length(uniThresholds));
%     
%     for i = fliplr(1:length(uniThresholds))
%         pointsBelow = zeros(size(freqSizeStrucutre));
%         
%         pointsBelow(freqSizeStrucutre <= uniThresholds(i) & freqSizeStrucutre > 0) = 1;
%         
%         indsBelow = find(pointsBelow);
%         
%         [xBelow, yBelow] = ind2sub(size(freqSizeStrucutre),  indsBelow);
%         
%         for j = 1:size(freqSizeStrucutre,2)
%             yInds = find(yBelow == j);
%             
%             if ~isempty(yInds)
%                 otherInds = find(yBelow ~= j);
% 
%                 for k = 1:length(yInds)
%                     c(j,i) = c(j,i) + length(find(xBelow(yInds(k)) == xBelow(otherInds)));
%                 end
%             else
%                 c(j,i) = NaN;
%             end
%         end
%     end
    
    % ceq is equality constraint, solves for ceq(x) = 0
    % calculate ceq for each size
    % find number of troughs in given size and aim for 1
    
    ceq = zeros(size(freqSizeStrucutre,2),1);
    
    allMinima = zeros(size(freqSizeStrucutre));
    
    for i = 1:size(freqSizeStrucutre,2)
        inactLevels = freqSizeStrucutre(:,i);
        
        inactInds = find(inactLevels > 0);
        
        inactLevels = inactLevels(inactInds);
        
        %%% Allow option to just count minima, not require points (broad peaks)...
        
        % Check not all equal
        if length(inactLevels) > sum(inactLevels == min(inactLevels))
            % Only gets interior minima, and not edges...
            % Calculate both flat and first, flat always needed to relate between sizes
            minimaPointsAll = islocalmin(inactLevels, 'FlatSelection', 'all');
            
            minimaPointsFirst = islocalmin(inactLevels, 'FlatSelection', 'first');
            
            % Get first edge
            for j = 1:(length(minimaPointsAll)-1)
                if inactLevels(j) < inactLevels(j+1)
                    minimaPointsAll(1:j) = 1;
                    
                    minimaPointsFirst(1) = 1;
                    break;
                elseif inactLevels(j) > inactLevels(j+1)
                   break; 
                end

                % if equal, will continue
            end

            % Get last edge
            for j = fliplr(2:(length(minimaPointsAll)))
                if inactLevels(j-1) > inactLevels(j)
                    minimaPointsAll(j:end) = 1;
                    
                    minimaPointsFirst(j) = 1;
                    break;
                elseif inactLevels(j-1) < inactLevels(j)
                   break; 
                end

                % if equal, will continue
            end

            % allows single global minima, but multiple local minima
    %         minimaPoints = find(inactLevels == min(inactLevels));
    
        else
            minimaPointsAll = ones(length(inactLevels),1);
            
            minimaPointsFirst = zeros(length(inactLevels),1);
            minimaPointsFirst(1) = 1;
            
        end
        
        allMinima(inactInds,i) = minimaPointsAll;
        
        % Subtract 1 so goal is 1 point
        % Type 0 is all (forces point), type 1 is first (allows broad)
        if minimaType
            ceq(i) = sum(minimaPointsFirst) - 1;
        else
            ceq(i) = sum(minimaPointsAll) - 1;
        end
    end
    
    % 2nd equality constraint - minima should be higher freq than former
    % and lower freq the latter
    
    ceq2 = zeros(size(freqSizeStrucutre,2),1);
    
    % This assumes that inactivation is rising from left to right
    
    for i = 1:size(freqSizeStrucutre,2)
        
        minimaInds = find(allMinima(:,i));
        
        %%% Need to check this constratin works well...
        
        if i > 1
           for j = 1:i-1
               formerInds = find(allMinima(:,j));

               if isempty( intersect(minimaInds, formerInds))
                   indDif = formerInds(1) - minimaInds;
                   % Former inds should be lower, so diff greater than zero
                   if all(indDif < 0)
                        ceq2(i) = ceq2(i) + sum(indDif < 0);
                   end
               end
           end
        end
        
        if i < size(freqSizeStrucutre,2)
            for j = i+1:size(freqSizeStrucutre,2)
                latterInds = find(allMinima(:,j));

                if isempty( intersect(minimaInds, latterInds))
                    indDif = minimaInds(1) - latterInds;
                    %Latter inds should be higher, so diff greater than zero
                    if all(indDif < 0)
                        ceq2(i) = ceq2(i) + sum(indDif < 0);
                    end
                end
            end
        end
        
    end
    
    ceq = [ceq' ceq2']';
end

