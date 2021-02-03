function [absorbtion, extinction, power, amplitude] = calculatesphereabsorbtion(frequencyRange_rad, resonance_rad, ...
    mass, qualityFactor, chargeDifference, number, area, drive, mediumPermitivity)
    % from Yang...Sun 2016

    if nargin == 9
        if ~isempty(mediumPermitivity)
           error('Flux equaiton needs to be modified') 
        end
    end
    
    % Constants  - should be same as the main
    VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)
    
    LIGHT_SPEED = 299792458; % m/s
    
    absorbtion = zeros(length(frequencyRange_rad), 1);
    
    extinction = zeros(length(frequencyRange_rad), 1);
    
    power = zeros(length(frequencyRange_rad), 1);
    
    amplitude = zeros(length(frequencyRange_rad), 1);
    
    %%% Should just be vectorized rather than in loop
    
    for iFreq = 1:length(frequencyRange_rad)
        % Eqn 7
        amplitude(iFreq) = chargeDifference*drive./(mass*...
            sqrt((resonance_rad^2 - frequencyRange_rad(iFreq)^2)^2 + ...
            (resonance_rad*frequencyRange_rad(iFreq)/qualityFactor).^2));

        % Eqn 9
        power(iFreq) = resonance_rad * ...
            frequencyRange_rad(iFreq)^2 * mass * ...
            amplitude(iFreq)^2 / (2*qualityFactor);

        % Modify this eqn if medium permitivity is defined
        powerFlux = 0.5*VACCUM_PERMITIVITY*LIGHT_SPEED*drive^2;

        % Eqn 10
        extinction(iFreq) = power(iFreq)/powerFlux;

        % Eqn 13 - solution for absorbtion
        absorbtion(iFreq) = (1-exp(-extinction(iFreq)*...
            number/area));

    end

end

