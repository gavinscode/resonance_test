% clear; close all; clc

cd('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Resonance_test_git/Data')

dataAbs2008 = readmatrix('Influenza_abs_2008.csv');

dataAbs2009 = readmatrix('Influenza_abs_2009.csv');

dataAbs2016 = readmatrix('Influenza_abs_2016.csv');

dataInactFreq2016 = readmatrix('Influenza_inact_freq_2016.csv');

dataInactPower2016 = readmatrix('Influenza_inact_int_2016.csv');

dataInactPower2014 = readmatrix('Influenza_inact_int_2014.csv');

dataDiameter2017 = readmatrix('Diameter_dist_2017.csv');

dataMass1984 = readmatrix('Mass_dist_1984.csv');


dataInactPower2016(end,2) = 100;

% Plot freq
f1 = figure; hold on;

plot(dataAbs2008(:,1), dataAbs2008(:,2), 'k-.', 'linewidth',2)
plot(dataAbs2009(:,1), dataAbs2009(:,2), 'k:', 'linewidth',2)
plot(dataAbs2016(:,1), dataAbs2016(:,2), 'k', 'linewidth',2)

set(gca,'TickDir','out', 'LineWidth',1, 'FontSize',15);
set(gcf, 'color', 'w');

title('Influenza absorbtion', 'FontSize', 20)
xlabel('Frequency (GHz)', 'FontSize',20)
ylabel('Absorption (%)', 'FontSize',20)
legend('2008 - H1N1','2009 - ?', '2016 - H3N2', 'FontSize',15)

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
plot(dataInactPower2014(:,1), dataInactPower2014(:,2), 'm:o', 'linewidth',2, 'markersize', 8)
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
legend('Inactivation:2016 - H3N2', 'Inactivation:2014 - H3N2', 'Limit - local, occupational', 'Limit - local, public', ...
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

%% Test symmetry
if 1
    
    figure;
    subplot(2,2,1); hold on
    % Flip bottom of each
    [~, maxInd] = max(dataAbs2008(:,2));
    plot(dataAbs2008(maxInd:end,1)-dataAbs2008(maxInd,1), dataAbs2008(maxInd:end,2), 'k', 'linewidth',2)
    plot(-(dataAbs2008(1:maxInd,1)-dataAbs2008(maxInd,1)), dataAbs2008(1:maxInd,2), 'r', 'linewidth',2)
    
    xlim([0 10]); title('2008')
    
    subplot(2,2,2); hold on
    [~, maxInd] = max(dataAbs2009(:,2));
    plot(dataAbs2009(maxInd:end,1)-dataAbs2009(maxInd,1), dataAbs2009(maxInd:end,2), 'k', 'linewidth',2)
    plot(-(dataAbs2009(1:maxInd,1)-dataAbs2009(maxInd,1)), dataAbs2009(1:maxInd,2), 'r', 'linewidth',2)
    
    xlim([0 10]); title('2009')
    
    subplot(2,2,3); hold on
    [~, maxInd] = max(dataAbs2016(:,2));
    plot(dataAbs2016(maxInd:end,1)-dataAbs2016(maxInd,1), dataAbs2016(maxInd:end,2), 'k', 'linewidth',2)
    plot(-(dataAbs2016(1:maxInd,1)-dataAbs2016(maxInd,1)), dataAbs2016(1:maxInd,2), 'r', 'linewidth',2)
  
    xlim([0 10]); title('2016')
    
    subplot(2,2,4); hold on
    [~, maxInd] = max(dataInactFreq2016(:,2));
    plot(dataInactFreq2016(maxInd:end,1)-dataInactFreq2016(maxInd,1), dataInactFreq2016(maxInd:end,2), 'k', 'linewidth',2)
    plot(-(dataInactFreq2016(1:maxInd,1)-dataInactFreq2016(maxInd,1)), dataInactFreq2016(1:maxInd,2), 'r', 'linewidth',2)
    
    xlim([0 10]); title('2016 Inact')
end

%% Test power log
if 0
    figure;
    subplot(2,2,1); hold on
    plot(dataInactPower2016(:,1), dataInactPower2016(:,2), 'm-o', 'linewidth',2, 'markersize', 8)
    plot(dataInactPower2016(:,1), -dataInactPower2016(:,2)+100, 'm:o', 'linewidth',2, 'markersize', 8)
    % Rough UV curve from cornavirus paper
    plot([0.01 0.5 1 2]*400, [99 10 1 0.01], 'c-o', 'linewidth',2, 'markersize', 8)
    plot([0.01 0.5 1 2]*400, -[99 10 1 0.01]+100, 'c:o', 'linewidth',2, 'markersize', 8)
    
    subplot(2,2,2); hold on
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), 'm-o', 'linewidth',2, 'markersize', 8)
    plot(log10(dataInactPower2016(:,1)), -dataInactPower2016(:,2)+100, 'm:o', 'linewidth',2, 'markersize', 8)
    
    plot(log10([0.01 0.5 1 2]*400), [99 10 1 0.01], 'c-o', 'linewidth',2, 'markersize', 8)
    plot(log10([0.01 0.5 1 2]*400), -[99 10 1 0.01]+100, 'c:o', 'linewidth',2, 'markersize', 8)
    
    subplot(2,2,3); hold on
    plot(dataInactPower2016(:,1), log10(dataInactPower2016(:,2)), 'm-o', 'linewidth',2, 'markersize', 8)
    plot(dataInactPower2016(:,1), log10(-[dataInactPower2016(1:3,2)' 99]+100), 'm:o', 'linewidth',2, 'markersize', 8)
    
    plot([0.01 0.5 1 2]*400, log10([100 10 1 0.01]), 'c-o', 'linewidth',2, 'markersize', 8)
    plot([0.01 0.5 1 2]*400, log10(-[99 10 1 0.01]+100), 'c:o', 'linewidth',2, 'markersize', 8)
    
    subplot(2,2,4); hold on
    plot(log10(dataInactPower2016(:,1)), log10(dataInactPower2016(:,2)), 'm-o', 'linewidth',2, 'markersize', 8)
    plot(log10(dataInactPower2016(:,1)), log10(-[dataInactPower2016(1:3,2)' 99]+100), 'm:o', 'linewidth',2, 'markersize', 8)
    
    plot(log10([0.01 0.5 1 2]*400), log10([100 10 1 0.01]), 'c-o', 'linewidth',2, 'markersize', 8)
    plot(log10([0.01 0.5 1 2]*400), log10(-[99 10 1 0.01]+100), 'c:o', 'linewidth',2, 'markersize', 8)
end
