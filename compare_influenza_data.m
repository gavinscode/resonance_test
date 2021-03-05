clear; close all; clc

cd('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Resonance_test_git/Data')

dataAbs2009 = readmatrix('Influenza_abs_2009.csv');

dataAbs2016 = readmatrix('Influenza_abs_2016.csv');

dataInactFreq2016 = readmatrix('Influenza_inact_freq_2016.csv');

dataInactInt2016 = readmatrix('Influenza_inact_int_2016.csv');

% Plot freq
figure; hold on;

plot(dataAbs2009(:,1), dataAbs2009(:,2), 'k:')
plot(dataAbs2016(:,1), dataAbs2016(:,2), 'k')

title('Influenza absorbtion')
xlabel('Frequency (GHz)')
ylabel('Absorption (%)')
legend('2009', '2016')

% Plot inactivation as im
freqRange = 5:0.2:15;
intRange = 0:20:1000;

inactImage = zeros(length(freqRange), length(intRange));

interpInactFreq = interp1(dataInactFreq2016(:,1), dataInactFreq2016(:,2),...
    freqRange, 'linear', 0);

interpInactInt = interp1(dataInactInt2016(:,1), dataInactInt2016(:,2),...
    intRange, 'linear', 0);

% Place int inativation
[~, freqInd] = min(abs(freqRange - 8));
inactImage(freqInd, :) = interpInactInt/100;

intRef = zeros(length(dataInactInt2016(:,1)),1);
for i = 1:length(dataInactInt2016(:,1))
   [~, intRef(i)] = min(abs(intRange - dataInactInt2016(i,1))); 
end

%Place freq inactivation
[~, intInd] = min(abs(intRange - dataInactInt2016(end,1)));
inactImage(:, intInd) = interpInactFreq/100;

freqRef = zeros(length(dataInactFreq2016(:,1)),1);
for i = 1:length(dataInactFreq2016(:,1))
   [~, freqRef(i)] = min(abs(freqRange - dataInactFreq2016(i,1)));
end

safeRef = zeros(length(intRange),4);
for i = 1:length(intRange)
   %%% note, limit is specified as spatial peak value of PD in controlled environments, not public
   %%% Sun divides by 2 to compare to average power?
   safePower = (200*(freqRange(i)/3).^0.2)/2;
   
  [~, safeRef(i,1)] = min(abs(intRange - safePower));
  
  safePower = (18.56*(freqRange(i))^0.699)/2;
   
  [~, safeRef(i,2)] = min(abs(intRange - safePower));
  
  safePower = (50)/2;
   
  [~, safeRef(i,3)] = min(abs(intRange - safePower));
  
  safePower = (10)/2;
   
  [~, safeRef(i,4)] = min(abs(intRange - safePower));
end

figure;
subplot(1,3,3)
imshow(inactImage');
hold on

% Plot measurment points
plot(freqRef, intInd, 'rx')
plot(freqInd, intRef, 'mo')

% Plot safe level
plot(1:51, safeRef(:,1), 'b')
plot(1:51, safeRef(:,2), 'g')
plot(1:51, safeRef(:,3), 'b--')
plot(1:51, safeRef(:,4), 'g--')

axis on
set(gca,'YTick', 1:10:51, 'YTickLabel', intRange(1:10:51), ...
    'XTick', 1:10:51, 'XTickLabel', freqRange(1:10:51))

xlabel('Frequency (GHz)')
ylabel('Power density (W/m^2)')

% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];

colormap(cols)

hBar = colorbar;
set(hBar, 'Ticks', 0:0.2:1, 'TickLabels', 0:20:100)

% Plot linear as well
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), 'r-x')
plot(dataAbs2016(:,1), dataAbs2016(:,2)/max(dataAbs2016(:,2))*100, 'k')
xlabel('Frequency (GHz)')
ylabel('Inactivation/Absorption (%)')
xlim([5 15]); ylim([0 101])
axis square

legend('Inactivation x freq','Absorption')

subplot(1,3,2)
plot(dataInactInt2016(:,1), dataInactInt2016(:,2), 'm-o')
xlabel('Power density (W/m^2)')
ylabel('Inactivation (%)')
line([1 1]*(200*(8/3)^0.2)/2, [0 100], 'color', 'b')
line([1 1]*(18.56*(8)^0.699)/2, [0 100], 'color', 'g')
line([1 1]*(50)/2, [0 100], 'color', 'b', 'linestyle', '--')
line([1 1]*(10)/2, [0 100], 'color', 'g', 'linestyle', '--')
ylim([0 1000]); ylim([0 101])
axis square

legend('Inactivation x power', 'Safe power - peak, occupational (8 GHz)', 'Safe power - peak, public (8 GHz)', ...
    'Safe power - average, occupational (8 GHz)', 'Safe power - average, public (8 GHz)')

%% Experiment with distance
warning('Get contour lines instead of plotting points')

LIGHT_SPEED = 299792458; % m/s

VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

distanceRange = 0:0.1:10;
    distanceTick = 1:10:101;

% powerRange = 10:10:1000; % Actually density, taken at 0.1 m  
powerRange = 0:0.1:20; % Taken at source
    powerTick = 1:10:201;
    
%%% This is for 18 degree horn antenna used in study
    % Get parameters for others to try
gainDBI = 30;
directionality = 360/gainDBI %Gain is roughly directionality
antennaGain = 10^(gainDBI/10); %20 % Convert from dBi to numeric

fieldThresholdLow = 10; % optimistic - originally 50, Sun thinks 68 is low!
fieldThresholdExpected = 87; % Sun expected

dielectricConstant = 67.4; % Water 8Ghz, real

powerLimitOccupational = (200*(8/3)^0.2)/2; % at 8 GHz
powerLimitPublic = (18.56*(8)^0.699)/2; % at 8 GHz

powerImage = zeros(length(distanceRange), length(powerRange));

for i = 1:length(distanceRange)
    for j = 1:length(powerRange)
        % From https://www.ahsystems.com/EMC-formulas-equations/field-intensity-calculation.php
%         S = Gt*Pt/4/pi/R^2
        powerImage(i,j) = powerRange(j)*antennaGain/(4*pi*distanceRange(i)^2);

        % Guess based on inverse square law
%         powerImage(i,j) = powerRange(j)*distanceRange(1)^2/distanceRange(i)^2;
    end
end

fieldImageAir = sqrt(powerImage*2/(LIGHT_SPEED*VACCUM_PERMITIVITY));

%%% Replace with refaction calc later
fieldImageWater = sqrt(powerImage*2/(LIGHT_SPEED*VACCUM_PERMITIVITY*sqrt(dielectricConstant)));

fieldImageWaterDrop = sqrt(powerImage*2/(LIGHT_SPEED*VACCUM_PERMITIVITY))*(3/(dielectricConstant+2));

figure;
subplot(3,2,1)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, air
inds = find(fieldImageAir > fieldThresholdLow);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'r.')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Air - optimistic threshold')

subplot(3,2,3)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, water
inds = find(fieldImageWater > fieldThresholdLow);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'r.')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Water - optimistic threshold')

subplot(3,2,5)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, water
inds = find(fieldImageWaterDrop > fieldThresholdLow);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'r.')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Drop - optimistic threshold')

subplot(3,2,2)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, air
inds = find(fieldImageAir > fieldThresholdExpected);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'r.')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Air - pessimistic threshold')

subplot(3,2,4)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, air
inds = find(fieldImageWater > fieldThresholdExpected);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'ro')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Water - pessimistic threshold')

subplot(3,2,6)
imshow(log10(powerImage)/3); hold on
% Get inactivated area - low, air
inds = find(fieldImageWaterDrop > fieldThresholdExpected);
[powerX, powerY] = ind2sub(size(powerImage), inds);
plot(powerY, powerX, 'ro')

axis on
set(gca,'YTick', distanceTick, 'YTickLabel', distanceRange(distanceTick), ...
    'XTick', powerTick, 'XTickLabel', powerRange(powerTick))
ylabel('Distance (m)')
xlabel('Source power (W)')

title('Drop - pessimistic threshold')

% Indicate dangerous area - public
inds = find(powerImage > powerLimitPublic);
[powerX, powerY] = ind2sub(size(powerImage), inds);
subplot(3,2,1); plot(powerY, powerX, 'gx')
subplot(3,2,2); plot(powerY, powerX, 'gx')
subplot(3,2,3); plot(powerY, powerX, 'gx')
% subplot(3,2,4); plot(powerY, powerX, 'gx')
subplot(3,2,5); plot(powerY, powerX, 'gx')
subplot(3,2,6); plot(powerY, powerX, 'gx')

% Indicate dangerous area - occupational
inds = find(powerImage > powerLimitOccupational);
[powerX, powerY] = ind2sub(size(powerImage), inds);
subplot(3,2,1); plot(powerY, powerX, 'bx')
subplot(3,2,2); plot(powerY, powerX, 'bx')
subplot(3,2,3); plot(powerY, powerX, 'bx')
subplot(3,2,4); plot(powerY, powerX, 'bx')
subplot(3,2,5); plot(powerY, powerX, 'bx')
% subplot(3,2,6); plot(powerY, powerX, 'bx')

%% Quick experiment with looking at field inside water droplet

dielectricConstant = 67;

% If radius << smaller than wavelength
fieldReduction_1 = 1/(3/(dielectricConstant+2))

% If similar size, take intensity from Mie scattering
m = sqrt(dielectricConstant);

%%% Get proper radial average - how does intenal e-field vary
% << ~= >>

x = pi*5e-6/(LIGHT_SPEED/(8*10^9));
nj=5*round(2+x+4*x.^(1/3))+160;
eValues = sqrt(mie_esquare(m, x, nj));

% Integral as in Q could be better reference, but already looks to be between extremes
fieldReduction_2 = 1./[max(eValues) mean(eValues) min(eValues)]

% These two are similar, for large constant, differ a lot for smaller
% Geometric case - for normal incidence, will decrease with angle
warning('calc of n should use complex modulus')
fieldReduction_3 = 1/(1-((1-sqrt(dielectricConstant))/(1+sqrt(dielectricConstant)))^2)

% From keeping power constant
fieldReduction_4 = 1/sqrt(1/sqrt(dielectricConstant))