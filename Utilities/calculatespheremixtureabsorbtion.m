function [absorbtion, extinction] = calculatespheremixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
    mass, qualityFactor, chargeDifference, number, area, drive, mediumPermitivity)

    nResonance = length(resonance_rad);
    
    % Check array lengths match
        % Should be number and mass for each resonance
    if nResonance ~= length(number) | nResonance ~= length(mass)
        error('array lengths incorrect')
    end
    
    if length(chargeDifference) > 1
        if nResonance ~= length(chargeDifference)
            error('array lengths incorrect')
        end
        chargeVaries = 1;
    else
        chargeVaries = 0;
    end
    
    if length(qualityFactor) > 1
       error('Array QF not jet implemented') 
    end
    
    absorbtion = zeros(length(frequencyRange_rad),1);
    
    extinction = zeros(length(frequencyRange_rad),1);
    
    for jDiameter = 1:nResonance
                       
        if chargeVaries
            chargeToUse = chargeDifference(jDiameter);
        else
            chargeToUse = chargeDifference;
        end

        [tempAbs, tempEx] = calculatesphereabsorbtion(...
            frequencyRange_rad, resonance_rad(jDiameter), mass(jDiameter), qualityFactor, chargeToUse, ...
            number(jDiameter), area, drive, mediumPermitivity);   
        
        % absorbtions multiply together
        absorbtion = 1 - (1 - absorbtion) .* (1 - tempAbs);
        
        % Extincions add by number
        extinction = extinction + tempEx*number(jDiameter);
    end
    % Divide by total number
    extinction = extinction/sum(number);
    
    % previously took from absorbtion, identical
%     extinction = -area/sum(number).*log(1-absorbtion);
end

