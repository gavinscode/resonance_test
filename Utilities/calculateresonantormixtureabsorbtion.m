function [absorbtion, extinction] = calculateresonantormixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
        intensity, qualityFactor, number, area, drive, mediumPermitivity)

    % Currently takes multiple sizes with fixed q, but will adapt to take
    % varying q and Q
 
    VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)
    
    LIGHT_SPEED = 299792458; % m/s

    nResonance = length(resonance_rad);

    if ~isempty(mediumPermitivity)
       error('Flux equaiton needs to be modified') 
    end

    if nResonance ~= length(number)
        error('array lengths incorrect')
    end
    
    if length(intensity) > 1
        if nResonance ~= length(intensity)
            error('intensity array length incorrect')
        end
        
        intensityVaries = 1;
    else
        intensityVaries = 0;
    end
    
    absorbtion = zeros(length(frequencyRange_rad),1);
    
%     figure; hold on
    
    % Calculate resonator power
    for jSize = 1:nResonance
        if intensityVaries
           intensityToUse = intensity(jSize); 
        else
           intensityToUse = intensity;
        end
        
        power = calculateresonatorpower(frequencyRange_rad, resonance_rad(jSize), intensityToUse, qualityFactor);
        
        tempEx = power/(0.5*VACCUM_PERMITIVITY*LIGHT_SPEED*drive^2);
        
        tempAbs = (1-exp(-tempEx*number(jSize)/area));
        
%         plot(tempAbs);
        
        absorbtion = 1 - (1 - absorbtion) .* (1 - tempAbs);
    end
    
    %%% Area should probably cancel out
    extinction = -area/sum(number).*log(1-absorbtion);
end

