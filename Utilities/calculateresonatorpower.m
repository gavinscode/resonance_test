function [power] = calculateresonatorpower(frequencyRange_rad, resonance_rad, ...
    intensity, gamma)
    % From Matthew Schwartz, Lecture 2 driven oscillators
    
    power = zeros(length(frequencyRange_rad), 1);
    
%     test = zeros(length(frequencyRange_rad), 1);
    
    for iFreq = 1:length(frequencyRange_rad)
        % drive and mass incorperated into intensity - drive^2/mass
        power(iFreq) = (intensity^2/(2*gamma))*(gamma*frequencyRange_rad(iFreq))^2/...
            ((resonance_rad^2-frequencyRange_rad(iFreq)^2)^2 + (gamma*frequencyRange_rad(iFreq))^2);
        
        % Test Lorentzian - not really equal ???
%         gammaT = gamma/1.5;
%         test(iFreq) = intensity^2/(pi*gammaT)*(gammaT^2/((frequencyRange_rad(iFreq) - resonance_rad)^2 + gammaT^2));
        
    end
        
end

