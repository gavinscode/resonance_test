function [absorbtion, extinction] = calculateresonantormixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
        intensity, qualityFactor, number, area, mediumPermitivity)

    % Currently takes multiple sizes with fixed q, but will adapt to take
    % varying q and Q
 
    VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)
    
    LIGHT_SPEED = 299792458; % m/s

    nResonance = length(resonance_rad);

    if ~isempty(mediumPermitivity)
       error('Flux equaiton needs to be modified') 
    end

    %%% If just getting extinction, could just use large value (10^18)
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
    
    if length(qualityFactor) > 1
        if nResonance ~= length(qualityFactor)
            error('QF array length incorrect')
        end
        
        qfVaries = 1;
    else
        qfVaries = 0;
    end
    
    absorbtion = zeros(length(frequencyRange_rad),1);
    
    extinction = zeros(length(frequencyRange_rad),1);

    % Calculate resonator power
    for jSize = 1:nResonance
        if intensityVaries
           intensityToUse = intensity(jSize); 
        else
           intensityToUse = intensity;
        end
        
        if qfVaries
           qfToUse = qualityFactor(jSize); 
        else
           qfToUse = qualityFactor;
        end
        
        % Assumes drive is 1, but it cancels for extinction and absorbtion
        power = calculateresonatorpower(frequencyRange_rad, resonance_rad(jSize), intensityToUse, qfToUse);
        
        tempEx = power/(0.5*VACCUM_PERMITIVITY*LIGHT_SPEED);
        
        tempAbs = (1-exp(-tempEx*number(jSize)/area));
        
        % Multiply absorbtions
        absorbtion = 1 - (1 - absorbtion) .* (1 - tempAbs);
        
        % Weighted addition of extincitons
        extinction = extinction + tempEx*number(jSize);
    end
   
    extinction = extinction/sum(number);
    
    % Previously calculated extinction from combined absorbtion
%     extinction = -area/sum(number).*log(1-absorbtion);
end

