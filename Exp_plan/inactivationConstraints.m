function [c, ceq] = inactivationConstraints(weightsVector, freqSizeStrucutre, weightsIndexes)

    freqSizeStrucutre(weightsIndexes) = weightsVector;

    % change to linear power for convience
    freqSizeStrucutre = 10.^(freqSizeStrucutre);
    
    % c is inequality constaint, solves for c(x) < 0
    % calculate c for each size
    
    c = []; zeros(size(freqSizeStrucutre,2),1);
    
%     for i = 1:size(freqSizeStrucutre,2)
%         
%     end
    
    % ceq is equality constraint, solves for ceq(x) = 0
    % calculate ceq for each size
    % find number of troughs in given size and aim for 1
    
    ceq = zeros(size(freqSizeStrucutre,2),1);
    
    for i = 1:size(freqSizeStrucutre,2)
        inactLevels = freqSizeStrucutre(:,i);
        
        inactLevels(inactLevels == 0) = [];
        
        % Check not all equal
        if length(inactLevels) > sum(inactLevels == min(inactLevels))
            % Only gets interior minima, and not edges...
            minimaPoints = islocalmin(inactLevels, 'FlatSelection', 'all');

%             minimaInds = find(minimaPoints);
%             
%             % Get any extra minima points
%             if ~isempty(minimaInds)
%                 for j = minimaInds
% 
%                     ind = j+1;
% 
%                     while inactLevels(ind) == inactLevels(j)
% 
%                         minimaPoints(ind) = 1;
% 
%                         ind = ind + 1;
%                     end
% 
%                 end
%             end
            
            % Get first edge
            for j = 1:(length(minimaPoints)-1)
                if inactLevels(j) < inactLevels(j+1)
                    minimaPoints(1:j) = 1;
                    break;
                elseif inactLevels(j) > inactLevels(j+1)
                   break; 
                end

                % if equal, will continue
            end

            % Get last edge
            for j = fliplr(2:(length(minimaPoints)))
                if inactLevels(j-1) > inactLevels(j)
                    minimaPoints(j:end) = 1;
                    break;
                elseif inactLevels(j-1) < inactLevels(j)
                   break; 
                end

                % if equal, will continue
            end

            % allows single global minima, but multiple local minima
    %         minimaPoints = find(inactLevels == min(inactLevels));
    
        else
            minimaPoints = ones(length(inactLevels),1);
        end
        
        % Subtract 1 so goal is 1 point
        ceq(i) = sum(minimaPoints) - 1;
    end
end

