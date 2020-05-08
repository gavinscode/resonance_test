function [output] = eulersolversecondorder(inputSignal, parameterList, timeStep)
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
    
    for iStep = 2:sampleLength
        ddY = (inputSignal(iStep-1)-systemSpring*output(iStep-1)-...
        systemDamping*dY)/systemMass;
    
        output(iStep) = output(iStep-1) + dY*timeStep;

        dY = dY + ddY*timeStep;
        
    end
end

