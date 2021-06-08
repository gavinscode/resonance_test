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
        
        % Cover full range - break handled by other constraint
        %%% Covering encouraged multiple minima, removed
%         minimaInds = minimaInds(1):minimaInds(end);
        
        % looking back
        if i > 1
           for j = 1:i-1
               formerInds = find(allMinima(:,j));
%                formerInds = formerInds(1):formerInds(end);
               
               % Get inds that don't overlap
               minimaIndsOut = minimaInds(~ismember(minimaInds, formerInds));
               formerIndsOut = formerInds(~ismember(formerInds, minimaInds));
               
               if ~isempty(minimaIndsOut)
                   
                   % got through each ind out in current
                   for k = 1:length(minimaIndsOut)
                       % get closest ind in former
                       [~, nearInd] = min(abs(minimaIndsOut(k) - formerInds));
                       
                       % Try just doing if minima levels are equal
                       if freqSizeStrucutre(minimaIndsOut(k),i) == freqSizeStrucutre(formerInds(nearInd),j)
                           distDiff = minimaIndsOut(k) - formerInds(nearInd);

                           % Former inds should be higher, so ok diff lower than zero
                           if distDiff > 0
                               % if not add to error
                               %%% Not sure which error model is best...
    %                            ceq2(i) = 1;
                               ceq2(i) = ceq2(i) + 1;
    %                            ceq2(i) = ceq2(i) + abs(distDiff);
                           end
                       end
                   end
               end
               
               %%% I don't think this will double count
               if ~isempty(formerIndsOut)
                   
                   % got through each ind out in former
                   for k = 1:length(formerIndsOut)
                       % get closest ind in current
                       
                       [~, nearInd] = min(abs(formerIndsOut(k) - minimaInds));
                       
                       if freqSizeStrucutre(minimaInds(nearInd),i) == freqSizeStrucutre(formerIndsOut(k),j)
                           distDiff = minimaInds(nearInd) - formerIndsOut(k);

                           % Former inds should be higher, so ok diff lower than zero
                           if distDiff > 0
                               % if not add to error
                               ceq2(i) = ceq2(i) + 1;
                           end
                       end
                   end
               end
           end
        end
        
        if i < size(freqSizeStrucutre,2)
            for j = i+1:size(freqSizeStrucutre,2)
                latterInds = find(allMinima(:,j));
%                 latterInds = latterInds(1):latterInds(end);
                
                % Get inds that don't overlap
                minimaIndsOut = minimaInds(~ismember(minimaInds, latterInds));
                latterIndsOut = latterInds(~ismember(latterInds, minimaInds));
                
                if ~isempty(minimaIndsOut)
                   
                   % got through each ind out in current
                   for k = 1:length(minimaIndsOut)
                       % get closest ind in former
                       
                       [~, nearInd] = min(abs(minimaIndsOut(k) - latterInds));
                       
                       if freqSizeStrucutre(minimaIndsOut(k),i) == freqSizeStrucutre(latterInds(nearInd),j)
                           distDiff = minimaIndsOut(k) - latterInds(nearInd);

                           % Latter inds should be lower, so ok diff higher than zero
                           if distDiff < 0
                               % if not add to error
                               ceq2(i) = ceq2(i) + 1;
                           end
                       end
                   end
               end
                
                if ~isempty(latterIndsOut)
                   
                   % got through each ind out in former
                   for k = 1:length(latterIndsOut)
                       % get closest ind in current
                       
                       [~, nearInd] = min(abs(latterIndsOut(k) - minimaInds));
                       
                       if freqSizeStrucutre(minimaInds(nearInd),i) == freqSizeStrucutre(latterIndsOut(k),j)
                           distDiff = minimaInds(nearInd) - latterIndsOut(k);

                           % Latter inds should be lower, so ok diff higher than zero
                           if distDiff < 0
                               % if not add to error
                               ceq2(i) = ceq2(i) + 1;
                           end
                       end
                   end
               end
            end
        end
        
    end
    
    ceq = [ceq' ceq2']';
end

