clear; clc
close all

testSymolically = 0

warning('Murray calcs only valid at mR = m/2')

%%% Take away
% Sun paper and solution give equal simulated and analytical power
    % Stimulated amplitude is also equal, but analytical reduces differently
% Sun eqns can be used confidently

% Murray solution gives same max power, but at different freq. 
    % Thus amplitude and phase differ at some freqs

if testSymolically
    syms k m b q E t w
    
    mR = m/2;

    % Given defs in Sun virus
    w0 = sqrt(k/mR);
    Q = w0*mR/b;

    w = w0;
    
    % Solves using modified relation from Murray - only correct at mR = m/2
    wM = 2*sqrt(k/m);
    
    % Get terms
    A_SunPaper = q*E/(mR*sqrt((w0^2 - w^2)^2 + (w0*w/Q)^2));

    %Does not work symbolically at resonance (w = w0) denominator goes to 0
    theta_SunPaper = atan( w0*w/(Q*(w0^2 - w^2)));

    A_SunSolution = (k - mR*w^2)/((k - mR*w^2)^2 + (w*b)^2);

    B_SunSolution = (w*b)/((k - mR*w^2)^2 + (w*b)^2);

    A_MurraySolution = (2*k - m*wM^2/2)/((2*k - m*wM^2/2)^2 + (2*wM*b)^2);
    
    B_MurraySolution = (2*wM*b)/((2*k - m*wM^2/2)^2 + (2*wM*b)^2);
    
    % Reduce x equations
    %%% Need to get a clearer expression for amplitude
    xMotion_SunSolution = q*E*(A_SunSolution*cos(w*t) + B_SunSolution*sin(w*t))

    xMotion_SunPaper = A_SunPaper*cos(w*t - theta_SunPaper)

    xMotion_MurraySolution = q*E*(A_MurraySolution*cos(w*t) + B_MurraySolution*sin(w*t))
    
    % Reduce power equations
    powerCycle_SunSolution = ((q*E)^2*w/2) * B_SunSolution

    powerCycle_SunPaper = (w0*w^2*mR/(2*Q)) * (A_SunPaper)^2

    powerCycle_MurraySolution = ((q*E)^2*wM/2) * B_MurraySolution * 2

    % Only at resonance (w = w0)
    powerCycle_MurrayPaper = (q*E)^2/(2*b)
else
    % Parameters directly from Sun virus
    m = 161*1.6605*10^-21; % Convert from MDa, full/non-hydrated
    
    q = 1.16*10^7*1.602176634*10^-19;
    
    E = 1;
    
    w0 = 8.22*10^9; 
    
    Q = 1.95;

    mR = m/2; %m*0.1*0.9;% 
    
    k = w0^2*mR;
    
    b = w0*mR/Q;
    
    w = (1:1:15)*10^9;
    
    % Time up to one cycle
    t = 0:10^-10:(2*pi/w(1));

    % Get terms
    A_SunPaper = q*E./(mR*sqrt((w0^2 - w.^2).^2 + (w0*w/Q).^2));

    %Does not work symbolically at resonance (w = w0) denominator goes to 0
    theta_SunPaper = atan( w0*w./(Q*(w0^2 - w.^2)));

    A_SunSolution = (k - mR*w.^2)./((k - mR*w.^2).^2 + (w*b).^2);

    B_SunSolution = (w*b)./((k - mR*w.^2).^2 + (w*b).^2);
    
    A_MurraySolution = (2*k - m.*w.^2/2)./((2*k - m*w.^2/2).^2 + (2.*w*b).^2);
    
    B_MurraySolution = (2.*w*b)./((2*k - m.*w.^2/2).^2 + (2.*w*b).^2);
    
    % Plot paths
    %%% Need to get a clearer expression for amplitude
    xMotion_SunSolution = q*E*(A_SunSolution(1)*cos(w(1)*t) + B_SunSolution(1)*sin(w(1)*t));

    xMotion_SunPaper = A_SunPaper(1)*cos(w(1)*t - theta_SunPaper(1));

    xMotion_MurraySolution = q*E*(A_MurraySolution(1)*cos(w(1)*t) + B_MurraySolution(1)*sin(w(1)*t));
    
    figure;
    
    subplot(2,2,1); hold on
    plot(t*10^9, xMotion_SunSolution, 'g')
    plot(t*10^9, xMotion_SunPaper, 'bx')
    plot(t*10^9, xMotion_MurraySolution, 'r')
    title('Position')
    
    subplot(2,2,3); hold on
    plot(t*10^9, xMotion_SunPaper-xMotion_SunSolution, 'k')
    title('Position Difference')
    
    % Plot power
    powerCycle_SunSolution = ((q*E)^2.*w/2) .* B_SunSolution;

    powerCycle_SunPaper = (w0*w.^2*mR/(2*Q)) .* (A_SunPaper).^2;

    powerCycle_MurraySolution = ((q*E)^2.*w/2) .* B_MurraySolution * 2;

    % Only at resonance (w = w0)
    powerCycle_MurrayPaper = (q*E)^2/(2*b);
    
    subplot(2,2,2); hold on
    plot(w/10^9, powerCycle_SunSolution, 'g')
    plot(w/10^9, powerCycle_SunPaper, 'xb')
    plot(w/10^9, powerCycle_MurraySolution, 'r')
    plot(w/10^9, powerCycle_MurrayPaper*ones(length(w),1), 'm')
    title('power')
    legend('Sun-Solution','Sun-Paper', 'Murray-Solution', 'Murray-Paper')

    subplot(2,2,4); hold on
    plot(w/10^9, powerCycle_SunPaper-powerCycle_SunSolution, 'k')
    title('Power Difference')
end