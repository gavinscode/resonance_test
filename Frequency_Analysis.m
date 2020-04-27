% Display frequncy components in different types of waves
clear

LIGHT_SPEED = 299792458;

% NearIR wavelength windows
window1Wavelength_m = [650 950]*10^-9;

window2Wavelength_m = [1100 1350]*10^-9;

window3Wavelength_m = [1600 1870]*10^-9;

window4Wavelength_m = [2100 2300]*10^-9;

%Target resonance
%%% Pulsed is very sensitive to this, 
%%% may be difference if this shifts between odd or even modes.
targetFrequency_hz = 8.22*10^9*1.25;

targetWavelength_m = LIGHT_SPEED/targetFrequency_hz;

targetPeriod_s = 1/targetFrequency_hz;

periodsToSimulate = 30;

simulationDuration_s = periodsToSimulate*targetPeriod_s;

% Convert wavelengths to frequncies
window1Frequency_hz = LIGHT_SPEED./window1Wavelength_m;

window2Frequency_hz = LIGHT_SPEED./window2Wavelength_m;

window3Frequency_hz = LIGHT_SPEED./window3Wavelength_m;

window4Frequency_hz = LIGHT_SPEED./window4Wavelength_m;

% Sampling rate
pointsPerCarrierWave = 10; %No qulatitative change 10 - 100

samplingFrequency_hz = window1Frequency_hz(1)*pointsPerCarrierWave;

samplingPeriod_s = 1/samplingFrequency_hz;

sampleLength = round(simulationDuration_s/samplingPeriod_s);



% Create signals
time_s = (0:sampleLength-1)*samplingPeriod_s; 

% Target signal is just sine wave
targetSignal = sin(2*pi*targetFrequency_hz*time_s);

% Ideal pulsed signal is half-rectified square wave
pulsedSignal = targetSignal;

pulsedSignal(pulsedSignal < 0) = 0;

pulsedSignal(pulsedSignal > 0) = 1;

% Add lightwave to pulses
carrierFrequency_hz = window4Frequency_hz(1);

lightPulseSignal = sin(2*pi*carrierFrequency_hz*time_s).*pulsedSignal;

% Create beating signal - subtract so longer wavelength
upperFrequncy_hz = carrierFrequency_hz - targetFrequency_hz;

beatingSignal = sin(2*pi*carrierFrequency_hz*time_s) + sin(2*pi*upperFrequncy_hz*time_s);

wavelengthDifference_nm = (LIGHT_SPEED/upperFrequncy_hz - LIGHT_SPEED/carrierFrequency_hz)*10^9;

% Create femto-pulsed signal
pulseDuration = 100*10^-15;

femtoPulsedSignal = pulsedSignal;

changeInds = find(diff(pulsedSignal) > 0) + 1;

% figure; hold on;
% plot(time_s, pulsedSignal) 
% plot(time_s(changeInds), pulsedSignal(changeInds),'rx') 

for i = 1:length(changeInds)-1
    inds = find(time_s > (time_s(changeInds(i)) + pulseDuration) & ... 
        time_s < time_s(changeInds(i+1)));
    
    femtoPulsedSignal(inds) = 0;
end

inds = find(time_s > (time_s(changeInds(end)) + pulseDuration));

femtoPulsedSignal(inds) = 0;

% Add lightwave to femto-pulses

lightFemtoPulsedSignal = sin(2*pi*carrierFrequency_hz*time_s).*femtoPulsedSignal;
    


% Frequency analysis 
frequencyReference = samplingFrequency_hz*(0:(sampleLength/2))/sampleLength;

% Target signal
targetPower = abs(fft(targetSignal)/sampleLength);

targetPower = targetPower(1:sampleLength/2+1);

targetPower(2:end-1) = 2*targetPower(2:end-1);

% Ideal pulsed
pulsedPower = abs(fft(pulsedSignal)/sampleLength);

pulsedPower = pulsedPower(1:sampleLength/2+1);

pulsedPower(2:end-1) = 2*pulsedPower(2:end-1);

% Lightwave pulsed
lightPulsedPower = abs(fft(lightPulseSignal)/sampleLength);

lightPulsedPower = lightPulsedPower(1:sampleLength/2+1);

lightPulsedPower(2:end-1) = 2*lightPulsedPower(2:end-1);

% Beating 
beatingPower = abs(fft(beatingSignal)/sampleLength);

beatingPower = beatingPower(1:sampleLength/2+1);

beatingPower(2:end-1) = 2*beatingPower(2:end-1);

% Femto-pulsed
femtoPulsedPower = abs(fft(femtoPulsedSignal)/sampleLength);

femtoPulsedPower = femtoPulsedPower(1:sampleLength/2+1);

femtoPulsedPower(2:end-1) = 2*femtoPulsedPower(2:end-1);

% Lightwave femto-pulsed
ligthFemtoPulsedPower = abs(fft(lightFemtoPulsedSignal)/sampleLength);

ligthFemtoPulsedPower = ligthFemtoPulsedPower(1:sampleLength/2+1);

ligthFemtoPulsedPower(2:end-1) = 2*ligthFemtoPulsedPower(2:end-1);



% Display results
%close all

figure; subplot(2,4,1); hold on

plot(time_s, targetSignal)

plot(time_s, pulsedSignal)

plot(time_s, femtoPulsedSignal)

subplot(2,4,2); hold on

plot(time_s, lightPulseSignal)

subplot(2,4,3); hold on

plot(time_s, beatingSignal)

title(sprintf('Wawvelength difference: %.f nm', wavelengthDifference_nm));

subplot(2,4,4); hold on

plot(time_s, lightFemtoPulsedSignal)

% Plot power spectra
subplot(2,4,5); hold on

plot(frequencyReference, targetPower)

plot(frequencyReference, pulsedPower)

plot(frequencyReference, femtoPulsedPower)

xlim([0 50*10^9])

subplot(2,4,6); hold on

plot(frequencyReference, lightPulsedPower)

xlim([0 50*10^9])

subplot(2,4,7); hold on

plot(frequencyReference, beatingPower)

xlim([0 50*10^9])

subplot(2,4,8); hold on

plot(frequencyReference, ligthFemtoPulsedPower)

xlim([0 50*10^9])