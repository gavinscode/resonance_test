function [output] = rksolversecondorder(inputSignal, parameterList, timeStep)
    % Note that strenght of input should be set in first argument
    
    % Assumes second order forced ODE of the form:
        % F = k*y + d*y' + m*y''
    
    if length(parameterList) ~= 3
       error('Incorrect number of paramteres') 
    end
    
    systemSpring = parameterList(1);
    
    systemDamping = parameterList(2);
    
    systemMass = parameterList(3);
    
    sampleLength = length(inputSignal);
    
    output = zeros(sampleLength, 1);
    
    dY = 0;
    
    % http://lampx.tugraz.at/~hadley/num/ch8/rk4ode2.php
    
    for iStep = 2:sampleLength
        
        error('Check integration vs Rung-Kutt-Nystrom in Kryszig');
        
        k1_1 = (inputSignal(iStep-1)-systemSpring*output(iStep-1)-...
            systemDamping*dY)/systemMass;
        
        k1_2 = dY;
        
        % Just average between adjacent time steps to get +dT/2 value
        
        k2_1 = ((inputSignal(iStep-1)+inputSignal(iStep))/2-...
            systemSpring*(output(iStep-1)+timeStep*k1_2/2)-...
            systemDamping*(dY+timeStep*k1_1/2))/systemMass;
        
        k2_2 = dY + timeStep*k1_2/2; 
        
        k3_1 = ((inputSignal(iStep-1)+inputSignal(iStep))/2-...
            systemSpring*(output(iStep-1)+timeStep*k2_2/2)-...
            systemDamping*(dY+timeStep*k2_1/2))/systemMass;
        
        k3_2 = dY + timeStep*k2_2/2; 
        
        k4_1 = (inputSignal(iStep)-...
            systemSpring*(output(iStep-1)+timeStep*k3_2)-...
            systemDamping*(dY+timeStep*k3_1))/systemMass;
        
        k4_2 = dY + timeStep*k3_2; 
        
        dY = dY + timeStep*(k1_1 + 2*(k2_1+k3_1) + k4_1)/6; 
        
        output(iStep) = output(iStep-1) + timeStep*(k1_2 + 2*(k2_2+k3_2) + k4_2)/6; 
        
    end
end

% k1_1 = (inputSignal(iStep-1)-systemSpring*output(iStep-1)-...
%             systemDamping*dY)/systemMass;
%         
% k2_1 = (mean(inputSignal(iStep-1:iStep))-systemSpring*output(iStep-1)-...
%     systemDamping*(dY+timeStep*k1_1/2))/systemMass;
% 
% k3_1 = (mean(inputSignal(iStep-1:iStep))-systemSpring*output(iStep-1)-...
%     systemDamping*(dY+timeStep*k2_1/2))/systemMass;
% 
% k4_1 = (inputSignal(iStep)-systemSpring*output(iStep-1)-...
%     systemDamping*(dY+timeStep*k3_1))/systemMass;
% 
% k1_2 = dY;
% 
% k2_2 = 
% 
% k3_2 = 
% 
% k4_2 =
% 
% output(iStep) = output(iStep-1) + timeStep*(k1_2 + 2*(k2_2+k3_2) + k4_2)/6; 
% 
% dY = dY + timeStep*(k1_1 + 2*(k2_1+k3_1) + k4_1)/6; 
