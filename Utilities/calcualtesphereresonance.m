function [frequency] = calcualtesphereresonance(radius, mode, l, n, soundSpeedL, soundSpeedT, guess)
    
    % Check sound is ratio
    if soundSpeedT < 1
        error('Ratio entered for tranverse speed?')
    end
    
    % n is indexed from zero.
    n = n + 1;
    
    if isempty(guess)
       guess = 10^10; 
    end
    
    % Only implemented for spherical functions with l = 1
    if mode == 'sph' & l == 1
        
        if n > 1 
           error('not yet implemented for higher order modes') 
        end
        
        angularFrequency = zeros(n,1);

        f = @(x)lambseqn_sph_1(x, radius, soundSpeedL, soundSpeedT);
        
        for iMode = 1:n
        
            %use optimizer to find fine intersection
            if iMode == 1
                % on first mode, just search from base
                angularFrequency(iMode) = fzero(f, guess*2*pi);
            else
                % for higher mode, search above last
                
                %%% Need a good strategy - maybe find each singularity with
                %%% fminbound and then search between them
                
%                 angularFrequency(iMode) = fzero(f, angularFrequency(iMode-1) + ...
%                     [1*10^9*2*pi angularFrequency(iMode-1)*4]);
            end
        end
    elseif mode == 'sph' & l ~= 1
        error('Only dipolar spherical mode is implemented')
    elseif mode ~= 'sph'
        error('No torsional modes are not implemented')
    end

    frequency = angularFrequency/2/pi;
    
end

% implements Eqn 1 from Sun 2015
function result = lambseqn_sph_1(w, r, cl, ct)

    if ~isempty(w)
        guess = w;
    end
    
    xi = (w*r/cl); %squiggle

    eta = (w*r/ct); %n
    
    %from wikipedia page on bessel functions
    j1_xi = sin(xi)/xi^2 - cos(xi)/xi;
    j2_xi = (3/xi^2 - 1)*sin(xi)/xi - 3*cos(xi)/xi^2;
    
    j1_eta = sin(eta)/eta^2 - cos(eta)/eta;
    j2_eta = (3/eta^2 - 1)*sin(eta)/eta - 3*cos(eta)/eta^2;
    
    % from Sun 2014
    result = 4*j2_xi/j1_xi*xi - eta^2 + ...
        2*j2_eta/j1_eta*eta;
end

