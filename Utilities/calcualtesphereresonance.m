function [frequency] = calcualtesphereresonance(radius, mode, l, n, soundSpeedL, soundSpeedT, startPoint, stepSize, plotFlag)
    % Guess should ideally be 1, but if first fmincon errors on a Nan,
    % it needs to be increased
    
    if isempty(startPoint)
       startPoint = 10^9*2*pi; 
    else
       startPoint = startPoint*2*pi; 
    end
    
    if isempty(stepSize)
       stepSize = 10^6*2*pi; 
    else
       stepSize = stepSize*2*pi; 
    end
    
    % Check sound is ratio
    if soundSpeedT < 1
        error('Ratio entered for tranverse speed?')
    end
    
    % n is indexed from zero.
    n = n + 1;
    
    % Plotting for debug
    if plotFlag
        freqz = (0:0.1:1000)*10^9;
        values = zeros(length(freqz),1);
        
        switch l
            case 0
                for i = 1:length(freqz)
                    values(i) = lambseqn_sph_0(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT)*200;
                end
            case 1
                for i = 1:length(freqz)
                    values(i) = lambseqn_sph_1(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT)*200;
                end
            case 2
                for i = 1:length(freqz)
                    values(i) = lambseqn_sph_2(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT);
                end
        end
        
        figure; 
        
        subplot(1,2,1); hold on;
        plot(freqz/10^9, values, 'g');
        ylim([-200 200])
        xlim([0 max(freqz)/10^9])
        
        subplot(1,2,2); hold on;
        plot(freqz/10^9, abs(values), 'g');
        ylim([0 200])
        xlim([0 max(freqz)/10^9])
    end
    
    % Only implemented for spherical functions, l = 0, 1, 2
    if mode == 'sph'
        % Select function appropriate for l
        switch l
            case 0
                % First function is normal
                    % Gets zero (crossings)
                f = @(x)(lambseqn_sph_0(x, radius, soundSpeedL, soundSpeedT));
                
                % Second function is  absolute value
                    % Get zeros (as minima)
                fabs = @(x)(abs(lambseqn_sph_0(x, radius, soundSpeedL, soundSpeedT)));
            case 1
                f = @(x)(lambseqn_sph_1(x, radius, soundSpeedL, soundSpeedT));
                
                fabs = @(x)(abs(lambseqn_sph_1(x, radius, soundSpeedL, soundSpeedT)));

            case 2
                f = @(x)(lambseqn_sph_2(x, radius, soundSpeedL, soundSpeedT));
                
                fabs = @(x)(abs(lambseqn_sph_2(x, radius, soundSpeedL, soundSpeedT)));
            otherwise
                error('Mode not implemented')
        end
        
        zeroFrequency = zeros(n,1);

        zerosFound = 0;
        
        oldValue = f(startPoint);
        
        its = 1;
        
        while zerosFound < n
            
            % Step up old value and find difference
            newValue = f(startPoint+stepSize);
            
            deriv = abs(newValue-oldValue);
            
            % Check if zero cross
                
            if ((oldValue < 0 & newValue > 0) | (oldValue > 0 & newValue < 0)) %& ...
                % And that difference is not large (indicates discontinuity)    
                if deriv < 50
             
                    % If so, find zero
                    zerosFound = zerosFound + 1;

                    % fminbnd seems on absolute seem to  work better than fzero on some points, unsure why...
                    % zeroFrequency(zerosFound) = fzero(f, startPoint);

                    zeroFrequency(zerosFound) = fminbnd(fabs, startPoint, startPoint+stepSize);

                    pointValue = f(zeroFrequency(zerosFound));

                    if plotFlag
                        %plot(zeroFrequency(zerosFound)/2/pi/10^9, -5, 'kx');

                        subplot(1,2,1);
                        plot(zeroFrequency(zerosFound)/2/pi/10^9, pointValue, 'ko');

                        subplot(1,2,2);
                        plot(zeroFrequency(zerosFound)/2/pi/10^9, pointValue, 'ko');
                    end
                    
                    its = 1;
                end
            end
            
            oldValue = newValue;
            
            startPoint = startPoint + stepSize;
         
            its = its + 1;
            
            if its > 10^10
               error('Not finding zeros, maybe problem with step size or derivative limit') 
            end
        end
    elseif mode ~= 'sph'
        error('No torsional modes are not implemented')
    end

    frequency = zeroFrequency/2/pi;
    
end

% implements Eqn 2 from Sun 2015 supplemental
function result = lambseqn_sph_0(w, r, cl, ct)
    % tested against Savoit app, r 4nm, 3700, 1540, Ok to 1GHz resolution

    xi = (w*r/cl); %greek squiggle

    %from wikipedia page on bessel functions (spherical bessel section)
    j0_xi = sin(xi)/xi;
    j1_xi = sin(xi)/xi^2 - cos(xi)/xi;
    
    result = 4*ct^2*j1_xi/(cl^2*xi)-j0_xi;
end

% implements Eqn 1 from Sun 2015
function result = lambseqn_sph_1(w, r, cl, ct)
    % tested against Savoit app, r 4nm, 3700, 1540, Ok to 1GHz resolution
    
    xi = (w*r/cl); %greek squiggle

    eta = (w*r/ct); %greek n
    
    %from wikipedia page on bessel functions (spherical bessel section)
    j1_xi = sin(xi)/xi^2 - cos(xi)/xi;
    j2_xi = (3/xi^2 - 1)*sin(xi)/xi - 3*cos(xi)/xi^2;
    
    j1_eta = sin(eta)/eta^2 - cos(eta)/eta;
    j2_eta = (3/eta^2 - 1)*sin(eta)/eta - 3*cos(eta)/eta^2;
    
    result = 4*j2_xi/j1_xi*xi - eta^2 + ...
        2*j2_eta/j1_eta*eta;
end

% implements Eqn 1 from Sun 2015 supplemental
function result = lambseqn_sph_2(w, r, cl, ct)
    % tested against Savoit app, r 4nm, 3700, 1540, Ok to 1GHz resolution

    xi = (w*r/cl); %greek squiggle

    eta = (w*r/ct); %greek n
    
    %from wikipedia page on bessel functions (spherical bessel section)
    j2_xi = (3/xi^2 - 1)*sin(xi)/xi - 3*cos(xi)/xi^2;
    j3_xi = (15/xi^3 - 6/xi)*sin(xi)/xi - (15/xi^2 - 1)*cos(xi)/xi;
    
    j2_eta = (3/eta^2 - 1)*sin(eta)/eta - 3*cos(eta)/eta^2;
    j3_eta = (15/eta^3 - 6/eta)*sin(eta)/eta - (15/eta^2 - 1)*cos(eta)/eta;
    
    l = 2;
    result = 4*(eta^2 + (l-1)*(l+2)*(eta*j3_eta/j2_eta-l-1))*xi*j3_xi/j2_xi - eta^4 + ...
        2*(l-1)*(2*l+1)*eta^2 + 2*(eta^2 - 2*l*(l-1)*(l+2))*eta*j3_eta/j2_eta;
end
