clear; close all; clc

cd('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Resonance_test_git/Data')

dataAbs2009 = readmatrix('Influenza_abs_2009.csv');

dataAbs2016 = readmatrix('Influenza_abs_2016.csv');

dataInactFreq2016 = readmatrix('Influenza_inact_freq_2016.csv');

dataInactPower2016 = readmatrix('Influenza_inact_int_2016.csv');

dataInactPower2016(end,2) = 100;

% Plot freq
f1 = figure; hold on;

plot(dataAbs2009(:,1), dataAbs2009(:,2), 'k:', 'linewidth',2)
plot(dataAbs2016(:,1), dataAbs2016(:,2), 'k', 'linewidth',2)

set(gca,'TickDir','out', 'LineWidth',1, 'FontSize',15);
set(gcf, 'color', 'w');

title('Influenza absorbtion', 'FontSize', 20)
xlabel('Frequency (GHz)', 'FontSize',20)
ylabel('Absorption (%)', 'FontSize',20)
legend('2009', '2016', 'FontSize',15)

% Plot inactivation as im
freqRange = 5:0.2:15;
powerRange = 0:20:1000;

inactImage = zeros(length(freqRange), length(powerRange));

interpInactFreq = interp1(dataInactFreq2016(:,1), dataInactFreq2016(:,2),...
    freqRange, 'linear', 0);

interpInactPower = interp1(dataInactPower2016(:,1), dataInactPower2016(:,2),...
    powerRange, 'linear', 0);

% Place power inativation
[~, freqInd] = min(abs(freqRange - 8));
inactImage(freqInd, :) = interpInactPower/100;

powerRef = zeros(length(dataInactPower2016(:,1)),1);
for i = 1:length(dataInactPower2016(:,1))
   [~, powerRef(i)] = min(abs(powerRange - dataInactPower2016(i,1))); 
end

%Place freq inactivation
[~, powerInd] = min(abs(powerRange - dataInactPower2016(end,1)));
inactImage(:, powerInd) = interpInactFreq/100;

freqRef = zeros(length(dataInactFreq2016(:,1)),1);
for i = 1:length(dataInactFreq2016(:,1))
   [~, freqRef(i)] = min(abs(freqRange - dataInactFreq2016(i,1)));
end

safeRef = zeros(length(powerRange),4);
for i = 1:length(powerRange)
   %%% note, limit is specified as spatial peak value of PD in controlled environments, not public
   %%% Sun divides by 2 to compare to average power?
   safePower = (200*(freqRange(i)/3).^0.2)/2;
   
  [~, safeRef(i,1)] = min(abs(powerRange - safePower));
  
  safePower = (18.56*(freqRange(i))^0.699)/2;
   
  [~, safeRef(i,2)] = min(abs(powerRange - safePower));
  
  safePower = (50)/2;
   
  [~, safeRef(i,3)] = min(abs(powerRange - safePower));
  
  safePower = (10)/2;
   
  [~, safeRef(i,4)] = min(abs(powerRange - safePower));
end

f2 = figure;
set(f2, 'WindowState', 'maximized')
pause(0.1)

subplot(1,4,3); hold on
imshow(inactImage');

% Plot measurment points
plot(freqInd, powerRef, 'mo', 'markersize', 8, 'linewidth', 2)
plot(freqRef, powerInd, 'rd', 'markersize', 8, 'linewidth', 2)
% Plot safe level
plot(1:51, safeRef(:,1), 'b', 'linewidth', 2)
plot(1:51, safeRef(:,2), 'g', 'linewidth', 2)
plot(1:51, safeRef(:,3), 'b--', 'linewidth', 2)
plot(1:51, safeRef(:,4), 'g--', 'linewidth', 2)

axis on
set(gca,'YTick', 1:10:51, 'YTickLabel', powerRange(1:10:51), ...
    'XTick', 1:10:51, 'XTickLabel', freqRange(1:10:51))
set(gca,'TickDir','out', 'LineWidth',1, 'FontSize',15);

title('Power by Frequency', 'FontSize', 20)
xlabel('Frequency (GHz)', 'FontSize', 20)
ylabel('Power density (W/m^2)', 'FontSize', 20)
xlim([1 51]); ylim([1 51])
% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];

colormap(cols)

% Plot linear as well
subplot(1,4,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), 'r-x', 'linewidth', 2, 'markersize', 8)
plot(dataAbs2016(:,1), dataAbs2016(:,2)/max(dataAbs2016(:,2))*100, 'k', 'linewidth', 1.5)
xlim([5 15]); ylim([0 100])
axis square 

set(gca,'TickDir','out', 'LineWidth',1, 'FontSize',15);

title('Inactivation by frequency', 'FontSize', 20)
xlabel('Frequency (GHz)', 'FontSize', 20)
ylabel('Inactivation/Absorption (%)', 'FontSize', 20)
legend('Inactivation','Absorption', 'FontSize', 20)

subplot(1,4,2); hold on
plot(dataInactPower2016(:,1), dataInactPower2016(:,2), 'm-o', 'linewidth',2, 'markersize', 8)
line([1 1]*(200*(8/3)^0.2)/2, [0 100], 'color', 'b', 'linewidth', 1.5)
line([1 1]*(18.56*(8)^0.699)/2, [0 100], 'color', 'g', 'linewidth', 1.5)
line([1 1]*(50)/2, [0 100], 'color', 'b', 'linestyle', '--', 'linewidth', 1.5)
line([1 1]*(10)/2, [0 100], 'color', 'g', 'linestyle', '--', 'linewidth', 1.5)
ylim([0 1000]); ylim([0 100])
axis square 

set(gca,'TickDir','out', 'LineWidth',1, 'FontSize',15);

title('Inactivation by power', 'FontSize', 20)
xlabel('Power density (W/m^2)', 'FontSize', 20)
ylabel('Inactivation (%)', 'FontSize', 20)
legend('Inactivation', 'Limit - local, occupational', 'Limit - local, public', ...
    'Limit - body, occupational', 'Limit - body, public', 'FontSize',15)

% Put color bar on seperate supbplot so it doesn't contract image
subplot(1,4,4); hold on
colormap(cols)
hBar = colorbar('west');
set(hBar, 'Ticks', 0:0.2:1, 'TickLabels', 0:20:100)
set(hBar,'TickDir','out', 'LineWidth',1, 'FontSize',15)
axis square
ylabel('Inactivation (%)', 'FontSize', 20)
axis off
temp = gca; temp.YLabel.Visible = 'on';
% temp = get(gca, 'Position');
% temp(1) = 0.77; set(gca, 'Position', temp)
% temp = get(hBar, 'Position');
% temp = [0.77 0.4, 0.008, 0.25]; set(hBar, 'Position', temp)

set(gcf, 'color', 'w');

%% Make simulated data

replicates = 3;

warning('Need to add random replicates with some variance')

% Phase 1.1 - also want to fit function to get peak 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 400;
simFreqTest_time = 1000;

curveMax = 80;
curveCenter = 8.5; 
curveSpread = 2.3;

% Phase 1.2 - indicate 5 points on funciton, indicate minimum effective power
simPowerTest_freqPoints = [0.01 0.5-(0.68*1.13)/2 0.5 0.5+(0.68*1.13)/2 0.99]; % from distribution
simPowerTest_power = 10:10:100;
simPowerTest_time = simFreqTest_time;

% test sigmoid and step-up
powerCentre = 70;
powerSlope = 0.08;
powerAmpOffset = [1 0.75 1 1 2]*0.6;

stepMin = 50*[1.4 1.2 1 0.9 0.8];

powerCols = cool(length(simPowerTest_freqPoints));

% Phase 1.3 - identify 50% and 100% effective durations
% Adjust times so it gives some more sense of variation
simTimeTest_freqPoints = simPowerTest_freqPoints; % from distribution
simTimeTest_effectPower = [1 2 4]; % multiplier
simTimeTest_times = [0.01 0.1 1 10 100]; % May start to get some band splitting issues.

timeSlope = 5;
timeSlopeOffset = [0.5 0.75 1.1 1.3 1.5];
powerSlopeOffset = [1 2 4];

% Phase 1.4
simScanTest_freqRange = simPowerTest_freqPoints([1 end]); % from distribution
simScanTest_freqSpacing = [1 2 4];
simScanTest_effectPower = [1 2 4]; % mult
simScanTest_minTime = [0.5 1]; % mult 

% Phase 1.5

% Phase 2.1

% Phase 2.2



% Make simulated plots
freqRange = min(simFreqTest_freqs):0.2:max(simFreqTest_freqs);
freqRangeFine = 0:0.05:max(simFreqTest_freqs);
powerRange = 0:5:500;

inactImage = zeros(length(freqRange), length(powerRange));

safeRef = zeros(length(freqRange),4);
for i = 1:length(freqRange)
  %%% note, limit is specified as spatial peak value of PD in controlled environments, not public
  safePower = (200*(freqRange(i)/3).^0.2)/2;
   
  [~, safeRef(i,1)] = min(abs(powerRange - safePower));
  
  safePower = (18.56*(freqRange(i))^0.699)/2;
   
  [~, safeRef(i,2)] = min(abs(powerRange - safePower));
  
  safePower = (50)/2;
   
  [~, safeRef(i,3)] = min(abs(powerRange - safePower));
  
  safePower = (10)/2;
   
  [~, safeRef(i,4)] = min(abs(powerRange - safePower));
end



% Phase 1.1 - frequency
%%% Should apply interpolated powerAmpOffset - but requires knowing limit freqs, so becomes circular...
simFreqTest_inact = curveMax*exp(-(simFreqTest_freqs-curveCenter).^2/(2*curveSpread^2));

interpInactFreq = interp1(simFreqTest_freqs, simFreqTest_inact,...
    freqRange, 'linear', 0);

interpInactFreq(interpInactFreq < 2) = 2;

%Place freq inactivation in map
[~, powerInd] = min(abs(powerRange - simFreqTest_power));
inactImage(:, powerInd) = interpInactFreq/100;

freqRef = zeros(length(simFreqTest_freqs),1);
for i = 1:length(simFreqTest_freqs)
   [~, freqRef(i)] = min(abs(freqRange - simFreqTest_freqs(i)));
end

figure;
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2)-(100-curveMax), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_freqs, simFreqTest_inact, 'rd', 'linewidth', 2, 'markersize', 8)
xlim([1 20]); ylim([0 100])

% fit curve
fittedCurve = fit(simFreqTest_freqs', simFreqTest_inact', 'gauss1');
curveInact = fittedCurve.a1*exp(-((freqRangeFine-fittedCurve.b1)/fittedCurve.c1).^2);

plot(freqRangeFine, curveInact, 'r-', 'linewidth', 2);

% take points for other tests
cdfCurve = cumsum(curveInact)/sum(curveInact);

% Power test
for i = 1:length(simPowerTest_freqPoints)
    inds = find(cdfCurve > simPowerTest_freqPoints(i));
    
    simPowerTest_freqPoints(i) = freqRangeFine(inds(1));
    
    plot(simPowerTest_freqPoints(i), fittedCurve.a1*exp(-((simPowerTest_freqPoints(i)-fittedCurve.b1)/fittedCurve.c1).^2), ...
        'o', 'markersize', 10, 'linewidth', 2, 'color', powerCols(i,:))
end

% Time test
for i = 1:length(simTimeTest_freqPoints)
    inds = find(cdfCurve > simTimeTest_freqPoints(i));
    
    simTimeTest_freqPoints(i) = freqRangeFine(inds(1));
end

% Scan test
for i = 1:length(simScanTest_freqRange)
    inds = find(cdfCurve > simScanTest_freqRange(i));
    
    simScanTest_freqRange(i) = freqRangeFine(inds(1));
end

% show image
subplot(1,3,2); hold on
imshow(inactImage');

plot(freqRef, powerInd, 'rd', 'markersize', 8, 'linewidth', 2)

% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];
colormap(cols)

plot(1:length(freqRange), safeRef(:,1), 'b', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,2), 'g', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,3), 'b--', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,4), 'g--', 'linewidth', 2)



% Phase 1.2 - power
figure; 

simPowerTest_inact_1 = zeros(length(simPowerTest_freqPoints), length(simPowerTest_power));
simPowerTest_inact_2 = zeros(length(simPowerTest_freqPoints), length(simPowerTest_power));

for i = 1:length(simPowerTest_freqPoints)
    % sigmoidal example
    curveVal = curveMax*exp(-(simPowerTest_freqPoints(i)-curveCenter).^2/(2*curveSpread^2)) * powerAmpOffset(i);
    
    simPowerTest_inact_1(i,:) = 1./(1+exp(-(simPowerTest_power-powerCentre)*powerSlope)) * curveVal;

    subplot(1,4,1); hold on
    plot(simPowerTest_power, simPowerTest_inact_1(i,:), 'o', 'color', powerCols(i,:), 'markersize', 6)
    
    % step-up example
    
    temp = simPowerTest_inact_1(i,:);
    
    temp(simPowerTest_power < stepMin(i)) = 0;
    
    simPowerTest_inact_2(i,:) = temp;
    
    subplot(1,4,2); hold on
    plot(simPowerTest_power, simPowerTest_inact_2(i,:), 'o', 'color', powerCols(i,:), 'markersize', 6)
end

subplot(1,4,1); hold on
plot(dataInactPower2016(:,1), dataInactPower2016(:,2), '-o', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([0 100]); ylim([0 50])

subplot(1,4,2); hold on
plot(dataInactPower2016(:,1), dataInactPower2016(:,2), '-', 'linewidth', 2, 'color', [0.5 0.5 0.5])
xlim([0 100]); ylim([0 50])

% Place sigmoidal in map
[~, minPowerInd] = min(abs(powerRange - simPowerTest_power(1)));
[~, maxPowerInd] = min(abs(powerRange - simPowerTest_power(end)));

freqIndArray = zeros(length(simPowerTest_freqPoints), 1);
powerRefArray = zeros(length(simPowerTest_power), 1);

for i = 1:length(simPowerTest_freqPoints)
    
    interpInactPower = interp1(simPowerTest_power, simPowerTest_inact_1(i,:),...
        powerRange(minPowerInd:maxPowerInd), 'linear', 0);

    interpInactPower(interpInactPower < 2) = 2;
    
    % Place power inativation
    [~, freqIndArray(i)] = min(abs(freqRange - simPowerTest_freqPoints(i)));
    inactImage(freqIndArray(i), minPowerInd:maxPowerInd) = interpInactPower/100;
    
end

for i = 1:length(simPowerTest_power)
   [~, powerRefArray(i)] = min(abs(powerRange - simPowerTest_power(i))); 
end

warning('Need to calculate minimum effective power given replicates')

% Temp solution
minEffPower = stepMin;

% Get index for maps
powerInd_MaxEff = zeros(length(simPowerTest_freqPoints), 3);

for i = 1:length(simPowerTest_freqPoints)
    [~, powerInd_MaxEff(i,1)] = min(abs(powerRange - minEffPower(i)));
    [~, powerInd_MaxEff(i,2)] = min(abs(powerRange - minEffPower(i)*simTimeTest_effectPower(2)));
    [~, powerInd_MaxEff(i,3)] = min(abs(powerRange - minEffPower(i)*simTimeTest_effectPower(3)));
end

% plot min eff on graphs - assuming first case
for i = 1:length(simPowerTest_freqPoints)

    curveVal = curveMax*exp(-(simPowerTest_freqPoints(i)-curveCenter).^2/(2*curveSpread^2)) * powerAmpOffset(i);
    
    temp = 1./(1+exp(-(minEffPower(i)-powerCentre)*powerSlope)) * curveVal;

    subplot(1,4,1); hold on
    plot(minEffPower(i), temp, '*', 'color', powerCols(i,:), 'markersize', 10)

end

% Show image
subplot(1,4,3); hold on
imshow(inactImage');

for i = 1:length(simPowerTest_freqPoints)
    for j = 1:length(simPowerTest_power)
        plot(freqIndArray(i), powerRefArray(j), 'o', 'markersize', 6, 'linewidth', 2, 'color', powerCols(i,:))
    end
end

plot(freqRef, powerInd, 'rd', 'markersize', 8, 'linewidth', 2)

for i = 1:length(simPowerTest_freqPoints)
    plot(freqIndArray(i), powerInd_MaxEff(i,1), '*', 'color', powerCols(i,:), 'markersize', 10)
end

% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];
colormap(cols)

plot(1:length(freqRange), safeRef(:,1), 'b', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,2), 'g', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,3), 'b--', 'linewidth', 2)
plot(1:length(freqRange), safeRef(:,4), 'g--', 'linewidth', 2)

% Plot finer image
subplot(1,4,4); hold on
imshow(inactImage');

for i = 1:length(simPowerTest_freqPoints)
    for j = 1:length(simPowerTest_power)
        plot(freqIndArray(i), powerRefArray(j), 'o', 'markersize', 6, 'linewidth', 2, 'color', powerCols(i,:))
    end
end

xlim([freqIndArray(1)-5 freqIndArray(end)+5])
ylim([powerRefArray(1)-5 max(powerInd_MaxEff(:))+5])

plot(freqIndArray, powerInd_MaxEff(:,1), 'color', [1 0.8 0], 'linewidth', 2)

for i = 1:length(simPowerTest_freqPoints)
    plot(freqIndArray(i), powerInd_MaxEff(i,1), '*', 'color', powerCols(i,:), 'markersize', 10, 'linewidth',2)
end
    
% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];
colormap(cols)



% Phase 1.3 - time
simTimesTest_inact_1 = zeros(length(simTimeTest_freqPoints), length(simTimeTest_effectPower), length(simTimeTest_times));

simTimesTest_inact_2 = zeros(length(simTimeTest_freqPoints), length(simTimeTest_effectPower), length(simTimeTest_times));

referenceInact = zeros(length(simTimeTest_freqPoints), length(simTimeTest_effectPower));

timeSizes = (4:4:(4*5));

% Plot absolute values
figAbsTime = figure;

for i = 1:length(simTimeTest_freqPoints)
    curveVal = curveMax*exp(-(simTimeTest_freqPoints(i)-curveCenter).^2/(2*curveSpread^2)) * powerAmpOffset(i);
    
    for j = 1:length(simTimeTest_effectPower)
        minPowerVal = 1./(1+exp(-(minEffPower(i)*simTimeTest_effectPower(j)-powerCentre)*powerSlope)) * curveVal;

        % Case 1, no effect of power
        simTimesTest_inact_1(i,j,:) = minPowerVal*(-1 + 2./(1+exp(-(simTimeTest_times)*timeSlope*timeSlopeOffset(i))));
        
        % Case 2, some effect of power
        simTimesTest_inact_2(i,j,:) = minPowerVal*(-1 + 2./(1+exp(-(simTimeTest_times)*timeSlope*timeSlopeOffset(i)*powerSlopeOffset(j))));
        
        referenceInact(i,j) = minPowerVal;
        
        for k = 1:length(simTimeTest_times)
            subplot(2,3,j); hold on
            plot(log10(simTimeTest_times(k)), permute(simTimesTest_inact_1(i,j,k),[3 1 2]), 's', 'color', powerCols(i,:), 'markersize', timeSizes(k))
            plot(log10(simPowerTest_time), minPowerVal, 'o', 'color', powerCols(i,:), 'markersize', 8);

            subplot(2,3, 3+j); hold on
            plot(log10(simTimeTest_times(k)), permute(simTimesTest_inact_2(i,j,k),[3 1 2]), 's', 'color', powerCols(i,:), 'markersize', timeSizes(k))
            plot(log10(simPowerTest_time), minPowerVal, 'o', 'color', powerCols(i,:), 'markersize', 8);
        end
    end
end

% Plot relative values
figRelTime = figure;
for i = 1:length(simTimeTest_freqPoints)
    for j = 1:length(simTimeTest_effectPower)
        for k = 1:length(simTimeTest_times)
            subplot(2,3,j); hold on
            plot(log10(simTimeTest_times(k)), permute(simTimesTest_inact_1(i,j,k),[3 1 2])/referenceInact(i,j), 's', 'color', powerCols(i,:), 'markersize', timeSizes(k))
            plot(log10(simPowerTest_time), 1, 'o', 'color', powerCols(i,:), 'markersize', 8);

            subplot(2,3, 3+j); hold on
            plot(log10(simTimeTest_times(k)), permute(simTimesTest_inact_2(i,j,k),[3 1 2])/referenceInact(i,j), 's', 'color', powerCols(i,:), 'markersize', timeSizes(k))
            plot(log10(simPowerTest_time), 1, 'o', 'color', powerCols(i,:), 'markersize', 8);
        end
    end
end

%%% Note, these are log plots, change x axis

%%% Get these from relative values - assume no effect of power for now
% Use points req for scanning;
effTime = zeros(length(simTimeTest_freqPoints), length(simScanTest_minTime));
effTime(:,1) = [7 5 3 1 0.5]; % 100%
effTime(:,2) = [7 5 3 1 0.5]*0.2; % 50%

timeMarkers = {'x', '+'};
timeLineWidth = [1 3];

% Sizes for min effective times - not really much change...
effTimeSizes = zeros(length(simTimeTest_times), length(simScanTest_minTime));
for i = 1:length(simScanTest_minTime)
    effTimeSizes(:,i) = interp1(log10(simTimeTest_times), timeSizes, log10(effTime(:,i)), 'linear', 0);
end

% Just plot on one figure for now
figure(figRelTime); subplot(2,3,1);
for i = 1:length(simTimeTest_freqPoints)
    for j = 1:length(simScanTest_minTime)
        plot(log10(effTime(i,j)), 1, timeMarkers{j}, 'color', powerCols(i,:), 'markersize', effTimeSizes(i,j), 'linewidth', timeLineWidth(j));
    end
end

% Represent on map with marker size

% Place inactivation on map for points not already covered 
%%% A bit sketch - not exactly the same format as other exp
%%% could also add to power sweep to fit curve better than linear (place on that figure)
for i = 1:length(simPowerTest_freqPoints)
    tempPowerInds = minPowerInd:powerInd_MaxEff(i,end);
    
    interpInactPower = interp1(minEffPower(i)*simTimeTest_effectPower, simTimesTest_inact_1(i,:,end),...
        powerRange(tempPowerInds), 'linear', 0);

    interpInactPower(interpInactPower < 2) = 2;
    
    % Place power inativation
    inds2Place = find(inactImage(tempPowerInds) == 0);
    
    inactImage(freqIndArray(i), tempPowerInds(inds2Place)) = interpInactPower(inds2Place)/100;
    
end

figure;
subplot(1,2,1); hold on
imshow(inactImage');

xlim([freqIndArray(1)-5 freqIndArray(end)+5])
ylim([powerRefArray(1)-5 max(powerInd_MaxEff(:))+5])

for i = 1:length(simTimeTest_freqPoints)
    for j = 1:length(simTimeTest_times)
        plot(freqIndArray(i), powerInd_MaxEff(i,1), 's', 'color', powerCols(i,:), 'markersize', timeSizes(j), 'linewidth',0.5)
        plot(freqIndArray(i), powerInd_MaxEff(i,2), 's', 'color', powerCols(i,:),  'markersize', timeSizes(j), 'linewidth', 0.5)
        plot(freqIndArray(i), powerInd_MaxEff(i,3), 's', 'color', powerCols(i,:), 'markersize', timeSizes(j), 'linewidth',0.5)
    end
end

for i = 1:length(simPowerTest_freqPoints)
    for j = 1:length(simPowerTest_power)
        plot(freqIndArray(i), powerRefArray(j), 'o', 'markersize', 6, 'linewidth', 0.5, 'color', powerCols(i,:))
    end
end

%%% May need better color ID here - link to 2D plots
plot(freqIndArray, powerInd_MaxEff(:,1), 'color', [1 0.8 0], 'linewidth', 2)
plot(freqIndArray, powerInd_MaxEff(:,2), 'color', [1 0.8 0], 'linewidth', 2)
plot(freqIndArray, powerInd_MaxEff(:,3), 'color', [1 0.8 0], 'linewidth', 2)

% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];
colormap(cols)

subplot(1,2,2); hold on
imshow(inactImage');

xlim([freqIndArray(1)-5 freqIndArray(end)+5])
ylim([powerRefArray(1)-5 max(powerInd_MaxEff(:))+5])

for i = 1:length(simTimeTest_freqPoints)
    % These look basically the same on size scaling
    for j = 1:length(simScanTest_minTime)
        plot(freqIndArray(i), powerInd_MaxEff(i,1), timeMarkers{j}, 'color', powerCols(i,:), 'markersize', effTimeSizes(i,j), 'linewidth',timeLineWidth(j))
        plot(freqIndArray(i), powerInd_MaxEff(i,2), timeMarkers{j}, 'color', powerCols(i,:), 'markersize', effTimeSizes(i,j), 'linewidth',timeLineWidth(j))
        plot(freqIndArray(i), powerInd_MaxEff(i,3), timeMarkers{j}, 'color', powerCols(i,:), 'markersize', effTimeSizes(i,j), 'linewidth',timeLineWidth(j))     
    end
end

for i = 1:length(simPowerTest_freqPoints)
    for j = 1:length(simPowerTest_power)
        plot(freqIndArray(i), powerRefArray(j), 'o', 'markersize', 6, 'linewidth', 0.5, 'color', powerCols(i,:))
    end
end

%%% May need better color ID here - link to 2D plots
plot(freqIndArray, powerInd_MaxEff(:,1), 'color', [1 0.8 0], 'linewidth', 2)
plot(freqIndArray, powerInd_MaxEff(:,2), 'color', [1 0.8 0], 'linewidth', 2)
plot(freqIndArray, powerInd_MaxEff(:,3), 'color', [1 0.8 0], 'linewidth', 2)

% Set up colors
cols = gray(101);
cols(1,:) = [0.2, 0, 0];
colormap(cols)



%% Phase 1.4 - scanning

inactivationCombos = zeros(length(simScanTest_freqSpacing), length(simScanTest_effectPower), length(simScanTest_minTime));

doseCombos = zeros(length(simScanTest_freqSpacing), length(simScanTest_effectPower), length(simScanTest_minTime));

durationCombos = zeros(length(simScanTest_freqSpacing), length(simScanTest_effectPower), length(simScanTest_minTime));

% do adjustment to bounds in place
adjustedFreqSpacing = zeros(length(simScanTest_freqSpacing),1);

for i = 1:length(simScanTest_freqSpacing)
    tempFreqs = simScanTest_freqRange(1):simScanTest_freqSpacing(i):simScanTest_freqRange(end);
    
    adjustedFreqSpacing(i) = (simScanTest_freqRange(end)-simScanTest_freqRange(1))/(length(tempFreqs));
    
    simScanTest_freqRange(1):adjustedFreqSpacing(i):simScanTest_freqRange(end)
end

% Also want to plot points tested on map...
for i = 1:length(simScanTest_freqSpacing)
    tempFreqs = simScanTest_freqRange(1):adjustedFreqSpacing(i):simScanTest_freqRange(end);
    
    tempPowerOffset = interp1(simPowerTest_freqPoints, powerAmpOffset, tempFreqs, 'linear', 0);
    
    curveVal = curveMax*exp(-(tempFreqs-curveCenter).^2/(2*curveSpread^2)) .* tempPowerOffset;
    
    % linear interpolate on powers
    tempPowers = interp1(simTimeTest_freqPoints, minEffPower, tempFreqs, 'linear', 0);
    
    for j = 1:length(simScanTest_effectPower)
        minPowerVal = 1./(1+exp(-(tempPowers*simScanTest_effectPower(j)-powerCentre)*powerSlope)) .* curveVal;
                    
        for k = 1:length(simScanTest_minTime) 
            %linear interp on times       
            %%% unfold this to full 2d when power effects time
%             timeInterpolant = scatteredInterpolant(simTimeTest_freqPoints',  minEffPower'*simTimeTest_effectPower(1), effTime(:,k), 'linear', 'none');
%             effTimeValues = timeInterpolant(tempFreqs, tempPowers*simScanTest_effectPower(j));
            
            effTimeValues = interp1(simTimeTest_freqPoints, effTime(:,k), tempFreqs, 'linear', 0);
            
            %%% No effect of power on time yet
            inactivationCombos(i,j,k) = sum(minPowerVal.*(-1 + 2./(1+exp(-(effTimeValues)*timeSlope*timeSlopeOffset(i)))));
    
            doseCombos(i,j,k) = wmean(tempPowers*simScanTest_effectPower(j), effTimeValues);
            
            durationCombos(i,j,k) = sum(effTimeValues);
        end
    end
end

inactivationCombos(inactivationCombos > 100) = 100;

powerEfficiencyCombos = inactivationCombos./doseCombos;

timeEfficiencyCombos = inactivationCombos./durationCombos;

figure;
% plot these results
for i = 1:length(simScanTest_effectPower)
   for j = 1:length(simScanTest_minTime) 

       subplot(3,length(simScanTest_minTime),j); hold on
       plot(adjustedFreqSpacing, inactivationCombos(:,i,j))

       subplot(3,length(simScanTest_minTime),j+2); hold on
       plot(adjustedFreqSpacing, doseCombos(:,i,j))
       
       subplot(3,length(simScanTest_minTime),j+4); hold on
       plot(adjustedFreqSpacing, durationCombos(:,i,j))
   end
end

%%% Also plot on map - just inactivation points?