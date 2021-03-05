function [power, amplitude] = calculateresonatorpower(frequencyRange_rad, resonance_rad, ...
    intensity, qualityFactor)
    % From Matthew Schwartz, Lecture 2 driven oscillators
    
    power = zeros(length(frequencyRange_rad), 1);
    
    amplitude = zeros(length(frequencyRange_rad), 1);
    
    % swapped to using consistent quality factor to match sphere resonator model
    gamma = resonance_rad/qualityFactor;
    
%     test = zeros(length(frequencyRange_rad), 1);
    
    %%% Should just be vectorized rather than in loop

    for iFreq = 1:length(frequencyRange_rad)
        % intensity currently represents intensity^2/mass
        % but will need to seperate when amplitude is implemented
        % eqn 34
        power(iFreq) = (intensity^2/(2*gamma))*(gamma*frequencyRange_rad(iFreq))^2/...
            ((resonance_rad^2-frequencyRange_rad(iFreq)^2)^2 + (gamma*frequencyRange_rad(iFreq))^2);
        
        % need to solve for specific solution for amplitude
        % roughly eqn 7 in Sun
        % will also need intesntiy and mass seperately  
        % eqn  24
%       amplitude(iFreq) = intensity/mass*()
        
        % Test Lorentzian - not really equal ???
%         gammaT = gamma/1.5;
%         test(iFreq) = intensity^2/(pi*gammaT)*(gammaT^2/((frequencyRange_rad(iFreq) - resonance_rad)^2 + gammaT^2));
        
    end
        
end

