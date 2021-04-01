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
