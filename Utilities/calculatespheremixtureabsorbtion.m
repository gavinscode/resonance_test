function [absorbtion, extinction] = calculatespheremixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
    mass, qualityFactor, chargeDifference, number, area, drive, mediumPermitivity)

    % Currently takes multiple sizes with fixed q, but will adapt to take
    % varying q and Q

    nResonance = length(resonance_rad);
    
    % Check array lengths match
        % Should be number and mass for each resonance
    if nResonance ~= length(number) | nResonance ~= length(mass)
        error('array lengths incorrect')
    end
    
    absorbtion = zeros(length(frequencyRange_rad),1);
    
    for jDiameter = 1:nResonance
                       
        tempAbs = calculatesphereabsorbtion(...
            frequencyRange_rad, resonance_rad(jDiameter), mass(jDiameter), qualityFactor, chargeDifference, ...
            number(jDiameter), area, drive, mediumPermitivity);    

        absorbtion = 1 - (1 - absorbtion) .* (1 - tempAbs);
    end
    
    % Predict extinction that would be recorded from measuring all spheres
    extinction = -area/sum(number).*log(1-absorbtion);

end

