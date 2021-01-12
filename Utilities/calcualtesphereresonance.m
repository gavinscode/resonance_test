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
    
    % To test
    figure; 
    freqz = (1:0.2:1000)*10^9;
    values = zeros(length(freqz),1);
    for i = 1:length(freqz)
        values(i) = lambseqn_sph_1(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT);
    end
    subplot(1,2,1); hold on;
    plot(freqz/10^9, values, 'b');
    subplot(1,2,2); hold on;
    plot(freqz/10^9, abs(values), 'b');
    
    values = zeros(length(freqz),1);
    for i = 1:length(freqz)
        values(i) = lambseqn_sph_0(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT);
    end
    subplot(1,2,1); hold on;
    plot(freqz/10^9, values*200, 'r');
    subplot(1,2,2); hold on;
    plot(freqz/10^9, abs(values)*200, 'r');
    
    values = zeros(length(freqz),1);
    for i = 1:length(freqz)
        values(i) = lambseqn_sph_2(freqz(i)*2*pi, radius, soundSpeedL, soundSpeedT);
    end
    subplot(1,2,1); hold on;
    plot(freqz/10^9, values, 'g');
    ylim([-200 200])
    subplot(1,2,2); hold on;
    plot(freqz/10^9, abs(values), 'g');
    ylim([0 200])
    
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
