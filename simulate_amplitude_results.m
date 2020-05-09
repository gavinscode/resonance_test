% From  Yang 2016, Figure 3B 
close all;
clear;
clc

%%% Note that results are ampltidue, but scaled to absorbtion

% Set up signals and time series
testFrequncies_hz = (5:0.1:15)*10^9;

testFrequncies_rad = testFrequncies_hz*2*pi;

testPeriods_s = 1./testFrequncies_hz;

periodsToSimulate = 10;

% Set duration based on longest period
simulationDuration_s = periodsToSimulate*max(testPeriods_s);

% Sampling rate
pointsPerWave = 100; 

samplingFrequency_hz = max(testFrequncies_hz)*pointsPerWave;

samplingPeriod_s = 1/samplingFrequency_hz;

sampleLength = round(simulationDuration_s/samplingPeriod_s);

time_s = (0:sampleLength-1)*samplingPeriod_s; 

% Measured paramteres - for reference if changed in model
measuredResonance_hz = 8.22*10^9;

measuredBandwidth = 4.22*10^9;

measuredQualityFactor = 1.95;

measuredAbsorbtion = 0.21;

reducedMass_kg = 14.5*1.6605*10^-21; % Convert from MDa to kg

% Model parameters
resonanceFrequency_hz = 8.3*10^9; 

resonantFrequency_rad = resonanceFrequency_hz*2*pi;

qualityFactor = 4.4;

systemSpring = resonantFrequency_rad^2*reducedMass_kg;

systemDamp = resonantFrequency_rad*reducedMass_kg/qualityFactor;

parameters = [systemSpring systemDamp reducedMass_kg];

% Simulate each using Runge-kutta
%%% Redo with FFT to get power as well

maxAmplitude = zeros(length(testFrequncies_hz),1);

figure; 
subplot(1,2,1); hold on
subplot(1,2,2); hold on

for iFreq = 1:length(testFrequncies_hz)
    % Create signal
    testSignal = makesinewave(testFrequncies_hz(iFreq), time_s);
    
    % Simualte virion motion
    outputMotion = rksolversecondorder(testSignal, parameters,...
        samplingPeriod_s);
    
    % Plot every 10th
    if rem(iFreq,10) == 0
        subplot(1,2,1); 
        plot(time_s, testSignal)

        subplot(1,2,2); 
        plot(time_s, outputMotion*10^12)
    end
    
    % Record result
    maxAmplitude(iFreq) = (max(outputMotion) - min(outputMotion))/2;
end

% Calculate analytical results
analyticalAmplitude = zeros(length(testFrequncies_hz),1);

%%% Amplitude does not match simulation perfectly - why not?

for iFreq = 1:length(testFrequncies_hz)
    % Eqn 7
    analyticalAmplitude(iFreq) = 1./(reducedMass_kg*...
        sqrt((resonantFrequency_rad^2 - testFrequncies_rad(iFreq).^2).^2 + ...
        (resonantFrequency_rad*testFrequncies_rad(iFreq)/qualityFactor).^2));

end

% Get peak 
[peakAmplitude, peakInd] = max(maxAmplitude);

simulatedPeakFrequency = testFrequncies_hz(peakInd)/10^9

% Scale amplitude
maxAmplitude = maxAmplitude/peakAmplitude*measuredAbsorbtion;

analyticalAmplitude = analyticalAmplitude/peakAmplitude*measuredAbsorbtion;

% Plot absorbtion spectrum
figure; hold on; ylim([0 0.25])
plot(testFrequncies_hz/10^9, maxAmplitude)

plot(testFrequncies_hz/10^9, analyticalAmplitude, '.-')

plot(testFrequncies_hz(peakInd)/10^9, maxAmplitude(peakInd), 'o')

% Indicate expected resonance
[~, resonanceInd] = min(abs(testFrequncies_hz - measuredResonance_hz));

plot(testFrequncies_hz(resonanceInd)/10^9, maxAmplitude(resonanceInd), 'x');

% Indicate bandwidth
upperBandInd = find(maxAmplitude(peakInd:end) < measuredAbsorbtion/2);

upperBandInd = upperBandInd(1) + peakInd - 1;

lowerBandInd = find(maxAmplitude(1:peakInd) < measuredAbsorbtion/2);

if ~isempty(lowerBandInd)
   lowerBandInd = lowerBandInd(end);
   
   plot(testFrequncies_hz([lowerBandInd upperBandInd])/10^9, maxAmplitude([lowerBandInd upperBandInd]), '-')
   
   simulatedBandwidth = (testFrequncies_hz(upperBandInd) - testFrequncies_hz(lowerBandInd))
   
   simulatedQualityFactor = simulatedPeakFrequency*10^9/simulatedBandwidth
else
    %Just use upper band as approximation
    plot(testFrequncies_hz(upperBandInd)/10^9, maxAmplitude(upperBandInd), '*')
    
    display('1 sided')
    
    simulatedBandwidth = (testFrequncies_hz(upperBandInd) - simulatedPeakFrequency*10^9)*2
   
    simulatedQualityFactor = simulatedPeakFrequency*10^9/simulatedBandwidth
end

ylabel('Amplitdue - scaled to absorbtion peak')
xlabel('Frequcny (GHz)')

