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
    
    for jDiameter = 1:nResonance
                       
        if chargeVaries
            chargeToUse = chargeDifference(jDiameter);
        else
            chargeToUse = chargeDifference;
        end

        tempAbs = calculatesphereabsorbtion(...
            frequencyRange_rad, resonance_rad(jDiameter), mass(jDiameter), qualityFactor, chargeToUse, ...
            number(jDiameter), area, drive, mediumPermitivity);   
        
        absorbtion = 1 - (1 - absorbtion) .* (1 - tempAbs);
    end
    
    % Predict extinction that would be recorded from measuring all spheres
    extinction = -area/sum(number).*log(1-absorbtion);

end

