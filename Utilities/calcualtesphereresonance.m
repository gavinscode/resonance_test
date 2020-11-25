function [frequency] = calcualtesphereresonance(radius, mode, l, n, soundSpeedL, speedRatio, guess)
    soundSpeedT = soundSpeedL*speedRatio;
    
    if isempty(guess)
       guess = 10^10; 
    end
     
%%% Check this
    %%% Need to get arbitrarily higher orders as well

%     testW = (0.1:0.1:200)*10^9*2*pi;
%     
%     temp = zeros(length(testW),1);
%     
%     for i = 1:length(testW)
%         temp(i) = lambseqn_sph_1(testW(i), radius, soundSpeedL, soundSpeedT);
%     end
%     
%     figure;
%     plot(testW/2/pi, temp); xlabel('Frequency (Hz)');
%     ylim([-20 20])
    
    % Only implemented for spherical functions with l = 1, n = 0
    if mode == 'sph' & l == 1
        
        %%% There are several, should search and pick lowest
        %n selects... (n zero is 1st)...
        
        %use optimizer to find fine intersection
        % f = @(x)lambseqn_sph_1(x, radius, soundSpeedL, soundSpeedT);
        %[angularFrequency] = fminsearch(f,guess*2*pi);

        angularFrequency = lambseqn_sph_1(guess*2*pi, radius, soundSpeedL, soundSpeedT);
        
    else
        error('mode not implemented')
    end

    frequency = angularFrequency/2/pi;
    
end

% implements Eqn 1 from Sun 2015
function result = lambseqn_sph_1(w, r, cl, ct)

    % Currently relies on guess a lot...
    % Better strategy could be to find inflection points, then search for zeros between them.

    if ~isempty(w)
        guess = w;
    end
    
    syms omega;
    
    % reduced frequncies
    xi = (omega*r/cl); %squiggle

    eta = (omega*r/ct); %n
    
%     xi = (w*r/cl); %squiggle
% 
%     eta = (w*r/ct); %n
    
    %from wikipedia page on bessel functions
    j1_xi = sin(xi)/xi^2 - cos(xi)/xi;
    j2_xi = (3/xi^2 - 1)*sin(xi)/xi - 3*cos(xi)/xi^2;
    
    j1_eta = sin(eta)/eta^2 - cos(eta)/eta;
    j2_eta = (3/eta^2 - 1)*sin(eta)/eta - 3*cos(eta)/eta^2;
    
    % from Sun 2014
%     result = 4*j2_xi/j1_xi*xi - eta^2 + ...
%         2*j2_eta/j1_eta*eta;
    
    eqn =  4*j2_xi/j1_xi*xi - eta^2 + ...
        2*j2_eta/j1_eta*eta == 0;
    
    result = vpasolve(eqn, omega, guess);
end

