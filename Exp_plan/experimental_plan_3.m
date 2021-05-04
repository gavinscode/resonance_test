compare_influenza_data
close all
%% - generate parameters for size - would be from SEM

influenzaSize_mean = 100*10^-9;
influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

influenzaVl = 1486; % Just matched in Saviot tool
influenzaVt = 743; %given x2 ratio

influenzaSize_samples = randn(500,1)*influenzaSize_std + influenzaSize_mean;

[influenzaSize_dist, influenzaSize_x] = hist(influenzaSize_samples,20);

influenzSize_resonances = zeros(length(influenzaSize_x),1);

for i = 1:length(influenzaSize_x)
    influenzSize_resonances(i) = calcualtesphereresonance(influenzaSize_x(i)/2, ...
                'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
end
 
%% third plan - sequential power steps.

nReps = 3;
absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% Phase 1.1 - also want to fit function to get peak 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 630;
simFreqTest_time = 1000;

curveMax = 100;
curveCenter = 8.5; 
curveSpread = 2.3;

% Phase 1.2
simPowerTest_inactLimit = 1; % from height of prediction bound, will be on min and max
% 1st freq points from inact distribution from 1st pass, 2nd from size dist for 2nd pass
    % Negative is below peak
    
%%% Adding setting so 0.75 taken from 1st curve and kept constant for 2nd step    
simPowerTest_freqPoints = {[1, 0.75, 0.95, -0.75, -0.95], []};
% not sure of best spacing to use here...
powerLogRange = 1:0.25:2.25;
simPowerTest_powers = round(10.^(powerLogRange))
simPowerTest_time = simFreqTest_time;

powerCols = cool(length(simPowerTest_powers));

% Just use single power threshold
% powerThreshold = 45;
% Use a interpolated threshold
powerThreshold = [15 -15]+45
powerThresholdFreqs = [-2.5 2.5] + 8.5;

useSharpStep = 1;
% need to adjust curve if not...

% Power decay
% Linear expression - note this is on log10 of power
powerLinearA = 0.709;
powerLinearB = -1.08;
zeroInactPoint = 10^(-powerLinearB/powerLinearA)


% Make simulated plots
% Set up image
freqRange = min(simFreqTest_freqs):0.2:max(simFreqTest_freqs);
freqRangeFine = 0:0.05:max(simFreqTest_freqs);
powerRange = round(10.^(0:0.05:3)); 0:10:1000;
powerRangeFine = round(10.^(0:0.005:3)); % 0:1:1000; 

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

%% Phase 1.1 - frequency
simFreqTest_inact = zeros(length(simFreqTest_freqs), nReps);
simFreqTest_inactRef = zeros(length(simFreqTest_freqs), nReps);
simFreqTest_FreqRef = zeros(length(simFreqTest_freqs), nReps);
simFreqTest_PowerRef = ones(length(simFreqTest_freqs), nReps)*simFreqTest_power;

% Get matched power across freq
if length(powerThreshold) == 1
    powerThresholdInterp_freqTest = powerThreshold*ones(length(simFreqTest_freqs),1);
else
    % Interp if more than one
    powerThresholdInterp_freqTest = interp1(powerThresholdFreqs, powerThreshold, simFreqTest_freqs, 'linear',0);
    powerThresholdInterp_freqTest(powerThresholdInterp_freqTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        simFreqTest_freqs(powerThresholdInterp_freqTest == 0), 'nearest','extrap');
end

for i = 1:length(simFreqTest_freqs)

    curveVal = curveMax*exp(-(simFreqTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    % Will always be above on this phase, but keep for consistancy
    if simFreqTest_power > powerThresholdInterp_freqTest(i)   
        if useSharpStep
            powerVal = powerLinearA*(log10(simFreqTest_power))+powerLinearB;
        else
            powerVal = powerLinearA*(log10(simFreqTest_power-powerThresholdInterp_freqTest(i)))+powerLinearB;
        end
        totalVal = powerVal * curveVal;
        totalVal(totalVal > 100) = 100;
        totalVal(totalVal < 0) = 0;
        
        simFreqTest_inactRef(i,:) = totalVal;
        simFreqTest_inact(i,:) = totalVal + absStd*randn(nReps,1);
    else
        simFreqTest_inact(i,:) = absStd*randn(nReps,1); %
    end

    simFreqTest_FreqRef(i,:) = simFreqTest_freqs(i);
end

% Can't be greater than 100
simFreqTest_inact(simFreqTest_inact > 100) = 100;

% fit curve
opts = fitoptions('gauss1', 'Lower', [0 1 0], 'Upper', [100 20 10]);

freqCurve = fit(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'gauss1', opts);

freq_curveInact = freqCurve(freqRangeFine);
freq_coefBound = confint(freqCurve);

freq_curveConfidence = predint(freqCurve, freqRangeFine, 0.95, 'Functional');
freq_curveObserved = predint(freqCurve, freqRangeFine, 0.95, 'obs');

% Get power test points
% From inactivation curve
inactCurveRatios = simPowerTest_freqPoints{1};
inactCurveFreqs = zeros(length(inactCurveRatios), 1);

for i = 1:length(inactCurveRatios)
   if inactCurveRatios(i) ~= 1
       freqInd = find(freq_curveInact/max(freq_curveInact) > abs(inactCurveRatios(i)));
       
       if inactCurveRatios(i) > 0
           inactCurveFreqs(i) = freqRangeFine(freqInd(end));
       else
           inactCurveFreqs(i) = freqRangeFine(freqInd(1));
       end
   else
       [~, freqInd] = max(freq_curveInact);
       
       inactCurveFreqs(i) = freqRangeFine(freqInd);
       
       centreFreq = freqRangeFine(freqInd);
   end
end

%%% Not using size curve now.
% Get freq curve from size with match spacing
% resonances_fine = interp1(influenzSize_resonances, influenzaSize_dist, freqRangeFine, 'linear', 0);
% resonances_fine(isnan(resonances_fine)) = 0;
% 
% sizeCurveRatios = simPowerTest_freqPoints{2};
% sizeCurveFreqs = zeros(length(sizeCurveRatios), 1);
% 
% for i = 1:length(sizeCurveRatios)
%    if sizeCurveRatios(i) ~= 1
%        freqInd = find(resonances_fine/max(resonances_fine) > abs(sizeCurveRatios(i)));
%        
%        if sizeCurveRatios(i) > 0
%            sizeCurveFreqs(i) = freqRangeFine(freqInd(end));
%        else
%            sizeCurveFreqs(i) = freqRangeFine(freqInd(1));
%        end
%    else
%        error('1 should not be in size points')
%    end
% end

simPowerTest_freqs = sort(inactCurveFreqs);

centrePowerInd = find(simPowerTest_freqs == centreFreq);

figure;
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'rd', 'linewidth', 2, 'markersize', 4)
xlim([1 20]); ylim([-10 110])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);
% plot(freqRangeFine, freq_curveObserved, 'r:', 'linewidth', 2)

% plot(simFreqTest_FreqRef(:,1), simFreqTest_inactRef(:,1), 'm-')

for i = 1:length(simPowerTest_freqs)
   line([1 1]*simPowerTest_freqs(i), [0 100], 'color', powerCols(i,:)); 
end

% show image
% Place freq inactivation in map
interpInactFreq = interp1(simFreqTest_freqs, mean(simFreqTest_inact,2),...
    freqRange, 'linear', 0);

interpInactFreq(interpInactFreq < 2) = 2;

[~, powerInd] = min(abs(powerRange - simFreqTest_power));
inactImage(:, powerInd) = interpInactFreq/100;

freqRef = zeros(length(simFreqTest_freqs),1);
for i = 1:length(simFreqTest_freqs)
   [~, freqRef(i)] = min(abs(freqRange - simFreqTest_freqs(i)));
end

% Show resonance distribution
plot(influenzSize_resonances, influenzaSize_dist, '-b', 'linewidth', 2)

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

%% Phase 1.2 - power

simPowerTest_inact = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_inactRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_freqRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_powerRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);

% Get matched roots across freq
if length(powerThreshold) == 1
    powerThresholdInterp_powerTest = powerThreshold*ones(length(simPowerTest_freqs),1);
else
    powerThresholdInterp_powerTest = interp1(powerThresholdFreqs, powerThreshold, simPowerTest_freqs, 'linear',0);
    powerThresholdInterp_powerTest(powerThresholdInterp_powerTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        simPowerTest_freqs(powerThresholdInterp_powerTest == 0), 'nearest','extrap');
end

for i = 1:length(simPowerTest_freqs)

    curveVal = curveMax*exp(-(simPowerTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simPowerTest_powers)

        if simPowerTest_powers(j) > powerThresholdInterp_powerTest(i)
            if useSharpStep
                powerVal = powerLinearA*(log10(simPowerTest_powers(j)))+powerLinearB;    
            else
                powerVal = powerLinearA*(log10(simPowerTest_powers(j)-powerThresholdInterp_powerTest(i)))+powerLinearB;    
            end
            
            totalVal = powerVal * curveVal;
            totalVal(totalVal > 100) = 100;     
            totalVal(totalVal < 0) = 0;

            simPowerTest_inactRef(i,j,:) = totalVal;
            
            simPowerTest_inact(i,j,:) = totalVal + absStd*randn(nReps,1);
        else
            simPowerTest_inact(i,j,:) = absStd*randn(nReps,1); %
        end
        
        simPowerTest_freqRef(i,j,:) = simPowerTest_freqs(i);
        simPowerTest_powerRef(i,j,:) = simPowerTest_powers(j);
    end
end

simPowerTest_inact(simPowerTest_inact > 100) = 100;

% Fit curves for each power level
figure; 

opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 0.5], 'Upper', [freqCurve.a1 freqCurve.b1 freqCurve.c1]);
ampLim = freqCurve.a1;
sigmaLim = freqCurve.c1;

amplitudeConfidence_1 = zeros(length(simPowerTest_powers),1);
powerCurveArray= cell(length(simPowerTest_powers),1);

for i = fliplr(1:length(simPowerTest_powers))
    
    subplot(1,4,1); hold on
    
    tempFreqs = simPowerTest_freqRef(:,i,:);
    tempInacts = simPowerTest_inact(:,i,:);
    
    plot(permute(tempFreqs, [1 3 2]), permute(tempInacts, [1 3 2]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols(i,:))

    if i ~= length(simPowerTest_powers)
        % Take limits from precedding
        
        opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 0.5], 'Upper', [ampLim freqCurve.b1 sigmaLim]);
    end
        
    powerCurve = fit(tempFreqs(:), tempInacts(:), 'gauss1', opts);
    power_coefBound = confint(powerCurve);

    powerCurveArray{i} = powerCurve;
    
    power_curveInact = powerCurve(freqRangeFine);
    
    plot(freqRangeFine, power_curveInact, 'color', powerCols(i,:))
        
    % plot height
    subplot(1,4,2); hold on
    plot(log10(simPowerTest_powers(i)), powerCurve.a1, 'o', 'color', powerCols(i,:));
    line([1 1]*log10(simPowerTest_powers(i)), [power_coefBound(1,1), power_coefBound(2,1)], 'color', powerCols(i,:))

    if ~isnan(power_coefBound(1,1))
        amplitudeConfidence_1(i) = power_coefBound(1,1);
        
        if power_coefBound(1,1) > 0
            ampLim = powerCurve.a1;
            
            sigmaLim = powerCurve.c1;
        end
    else
        amplitudeConfidence_1(i) = 0;
    end
    
    subplot(1,4,4); hold on
    plot(log10(simPowerTest_powers(i)), powerCurve.a1, 'o', 'color', powerCols(i,:));
    line([1 1]*log10(simPowerTest_powers(i)), [power_coefBound(1,1), power_coefBound(2,1)], 'color', powerCols(i,:))
    
    if ~isnan(power_coefBound(1,1))
        amplitudeConfidence_1(i) = power_coefBound(1,1);
    else
        amplitudeConfidence_1(i) = 0;
    end
        
%     subplot(3,7,4); hold on
%     plot(simPowerTest_freqRef(:,i,1), simPowerTest_inactRef(:,i,1), 'color', powerCols(i,:));
    
end

subplot(1,4,1); hold on

% plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
xlim([0 20]); ylim([-10 110])

subplot(1,4,2); hold on
for i = 1:length(powerThreshold)
    line([1 1]*log10(powerThreshold(i)), [0 50], 'color', 'k', 'linestyle', ':')
end
line([1 1]*log10(mean(powerThreshold)), [0 100], 'color', 'k', 'linestyle', ':')

line([1 1]*log10(zeroInactPoint), [0 100], 'color', 'm')

plot(log10(simFreqTest_power), freqCurve.a1, 'o', 'color', 'r');
line([1 1]*log10(simFreqTest_power), [freq_coefBound(1,1), freq_coefBound(2,1)], 'color', 'r')

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([0 3]); ylim([0 100])

% get points to test
% for power, first with above zero confidence bound
powerInd = find(amplitudeConfidence_1 > 1);

powerInd = powerInd(1);

simPowerTest_powers_2 = round(10.^(powerLogRange(powerInd-1):0.05:powerLogRange(powerInd)))

powerCols_2 = winter(length(simPowerTest_powers_2));

% for freq - move into 50% point from final curve
if ~isempty(simPowerTest_freqPoints{2})
    % Get curve again
    tempCurve = powerCurveArray{powerInd};
    power_curveInact = tempCurve(freqRangeFine);
    
    inactCurveRatios = simPowerTest_freqPoints{2};
    inactCurveFreqs = zeros(length(inactCurveRatios), 1);

    for i = 1:length(inactCurveRatios)
       if inactCurveRatios(i) ~= 1
           freqInd = find(power_curveInact/max(power_curveInact) > abs(inactCurveRatios(i)));

           if inactCurveRatios(i) > 0
               inactCurveFreqs(i) = freqRangeFine(freqInd(end));
           else
               inactCurveFreqs(i) = freqRangeFine(freqInd(1));
           end
       else
           [~, freqInd] = max(freq_curveInact);

           inactCurveFreqs(i) = freqRangeFine(freqInd);
       end
    end

    simPowerTest_freqs_2 = sort(inactCurveFreqs);
else
    % if empty, use the same
    simPowerTest_freqs_2 = simPowerTest_freqs;
end

% subplot(1,4,1);
% for i = 1:length(simPowerTest_freqs)
%    line([1 1]*simPowerTest_freqs(i), [0 100], 'color', powerCols(i,:)); 
% end
% 
% for i = 1:length(simPowerTest_freqs_2)
%    line([1 1]*simPowerTest_freqs_2(i), [0 100], 'color', powerCols_2(i,:)); 
% end

% plot(influenzSize_resonances, influenzaSize_dist, '-b', 'linewidth', 2)

%% Do second step of power test
%%% Set based on STD noise from 1st curve
nReps = 5;

simPowerTest_inact_2 = zeros(length(simPowerTest_freqs_2), length(simPowerTest_powers_2), nReps);
simPowerTest_inactRef_2 = zeros(length(simPowerTest_freqs_2), length(simPowerTest_powers_2), nReps);
simPowerTest_freqRef_2 = zeros(length(simPowerTest_freqs_2), length(simPowerTest_powers_2), nReps);
simPowerTest_powerRef_2 = zeros(length(simPowerTest_freqs_2), length(simPowerTest_powers_2), nReps);

% Get matched roots across freq
if length(powerThreshold) == 1
    powerThresholdInterp_powerTest = powerThreshold*ones(length(simPowerTest_freqs_2),1);
else
    powerThresholdInterp_powerTest = interp1(powerThresholdFreqs, powerThreshold, simPowerTest_freqs_2, 'linear',0);
    powerThresholdInterp_powerTest(powerThresholdInterp_powerTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        simPowerTest_freqs_2(powerThresholdInterp_powerTest == 0), 'nearest','extrap');
end

for i = 1:length(simPowerTest_freqs_2)

    curveVal = curveMax*exp(-(simPowerTest_freqs_2(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simPowerTest_powers_2)

        if simPowerTest_powers_2(j) > powerThresholdInterp_powerTest(i)
            if useSharpStep
                powerVal = powerLinearA*(log10(simPowerTest_powers_2(j)))+powerLinearB;    
            else
                powerVal = powerLinearA*(log10(simPowerTest_powers_2(j)-powerThresholdInterp_powerTest(i)))+powerLinearB;    
            end
            
            totalVal = powerVal * curveVal;
            totalVal(totalVal > 100) = 100;     
            totalVal(totalVal < 0) = 0;

            simPowerTest_inactRef(i,j,:) = totalVal;
            
            simPowerTest_inact_2(i,j,:) = totalVal + absStd*randn(nReps,1);
        else
            simPowerTest_inact_2(i,j,:) = absStd*randn(nReps,1); %
        end
        
        simPowerTest_freqRef_2(i,j,:) = simPowerTest_freqs_2(i);
        simPowerTest_powerRef_2(i,j,:) = simPowerTest_powers_2(j);
    end
end

simPowerTest_inact_2(simPowerTest_inact_2 > 100) = 100;

% Fit curves for each power level
% Use limit from one step up on previous
tempCurve = powerCurveArray{powerInd+1};
opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 0.5], 'Upper', [tempCurve.a1 freqCurve.b1 tempCurve.c1]);
ampLim = tempCurve.a1;
sigmaLim = tempCurve.c1;

amplitudeConfidence_2 = zeros(length(simPowerTest_powers_2),1);
powerCurveArray_2 = cell(length(simPowerTest_powers_2),1);

for i = fliplr(1:length(simPowerTest_powers_2))
    
    subplot(1,4,3); hold on
    
    tempFreqs = simPowerTest_freqRef_2(:,i,:);
    tempInacts = simPowerTest_inact_2(:,i,:);
    
    plot(permute(tempFreqs, [1 3 2]), permute(tempInacts, [1 3 2]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols_2(i,:))

    if i ~= length(simPowerTest_powers_2)
        % Take limits from precedding
        
        opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 0.5], 'Upper', [ampLim freqCurve.b1 sigmaLim]);
    end
        
    powerCurve = fit(tempFreqs(:), tempInacts(:), 'gauss1', opts);
    power_coefBound = confint(powerCurve);

    powerCurveArray_2{i} = powerCurve;
    
    power_curveInact = powerCurve(freqRangeFine);
    
    plot(freqRangeFine, power_curveInact, 'color', powerCols_2(i,:))

    % plot height
    subplot(1,4,4); hold on
    plot(log10(simPowerTest_powers_2(i)), powerCurve.a1, 'o', 'color', powerCols_2(i,:));
    line([1 1]*log10(simPowerTest_powers_2(i)), [power_coefBound(1,1), power_coefBound(2,1)], 'color', powerCols_2(i,:))
    
    if ~isnan(power_coefBound(1,1))
        amplitudeConfidence_2(i) = power_coefBound(1,1);
        
        if power_coefBound(1,1) > 0
            ampLim = powerCurve.a1;

            sigmaLim = powerCurve.c1;
        end
    else
        amplitudeConfidence_2(i) = 0;
    end

end

subplot(1,4,3); hold on

tempCurve = powerCurveArray{powerInd};
power_curveInact = tempCurve(freqRangeFine);
plot(freqRangeFine, power_curveInact, 'color', powerCols(powerInd,:))

plot(influenzSize_resonances, influenzaSize_dist/max(influenzaSize_dist)*tempCurve.a1, '-k', 'linewidth', 2)

% plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
% plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
xlim([5 8.5+3.5]); ylim([-10 30])

subplot(1,4,4); hold on
for i = 1:length(powerThreshold)
    line([1 1]*log10(powerThreshold(i)), [0 20], 'color', 'k', 'linestyle', ':')
end
line([1 1]*log10(mean(powerThreshold)), [0 40], 'color', 'k', 'linestyle', ':')

line([1 1]*log10(zeroInactPoint), [0 100], 'color', 'm')

plot(log10(simFreqTest_power), freqCurve.a1, 'o', 'color', 'r');
line([1 1]*log10(simFreqTest_power), [freq_coefBound(1,1), freq_coefBound(2,1)], 'color', 'r')

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 2]); ylim([0 40])

%% Add symmetry comparison
figure;

for i = 1:length(simPowerTest_freqs)
    subplot(2,5,i); hold on
    
    tempPowers = simPowerTest_powerRef(i,:,:);
    tempInacts = simPowerTest_inact(i,:,:);
    
    for j = 1:length(simPowerTest_powers)
        plot(log10(permute(tempPowers(1,j,:), [3 2 1])), permute(tempInacts(1,j,:), [3 2 1]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols(j,:))
    end
    
    if i == centrePowerInd
        plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    end
    
    plot(log10(permute(tempPowers(1,:,1), [3 2 1])), mean(permute(tempInacts(1,:,:), [3 2 1])), '-', 'linewidth', 2, 'markersize', 4, 'color', 'c')
    
    line([1 1]*log10(powerThresholdInterp_powerTest(i)), [0 100], 'color', 'k')
    
end

for i = 1:length(simPowerTest_freqs_2)
    subplot(2,5,i);
    
    tempPowers = simPowerTest_powerRef_2(i,:,:);
    tempInacts = simPowerTest_inact_2(i,:,:);
    
    for j = 1:length(simPowerTest_powers_2)
        plot(log10(permute(tempPowers(1,j,:), [3 2 1])), permute(tempInacts(1,j,:), [3 2 1]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols_2(j,:))
    end
    
    plot(log10(permute(tempPowers(1,:,1), [3 2 1])), mean(permute(tempInacts(1,:,:), [3 2 1])), '-', 'linewidth', 2, 'markersize', 4, 'color', 'g')
end

%%% Could combine duplicated powers in mean

for i = 1:centrePowerInd-1
    subplot(2,5,i+5); hold on
    
    % for 1st step
    topPowers = simPowerTest_powerRef(length(simPowerTest_freqs)-i+1,:,:);
    topInacts = simPowerTest_inact(length(simPowerTest_freqs)-i+1,:,:);
    
    bottomPowers = simPowerTest_powerRef(i,:,:);
    bottomInacts = simPowerTest_inact(i,:,:);
    
    meanInact = zeros(1,length(simPowerTest_powers));
    
    for j = 1:length(simPowerTest_powers)
        meanInact(j) = mean([permute(topInacts(1,j,:), [2 3 1]) permute(bottomInacts(1,j,:), [2 3 1])]);
        
        plot(log10(permute(topPowers(1,j,:), [3 2 1])), permute(topInacts(1,j,:), [3 2 1])-meanInact(j), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'r');
        
        plot(log10(permute(bottomPowers(1,j,:), [3 2 1])), permute(bottomInacts(1,j,:), [3 2 1])-meanInact(j), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'b');
    end
    
    plot(log10(permute(topPowers(1,:,1), [3 2 1])), mean(permute(topInacts(1,:,:), [3 2 1]))-meanInact, '-', 'linewidth', 2, 'markersize', 4, 'color', 'r')
    
    plot(log10(permute(bottomPowers(1,:,1), [3 2 1])), mean(permute(bottomInacts(1,:,:), [3 2 1]))-meanInact, '-', 'linewidth', 2, 'markersize', 4, 'color', 'b')
    
    % for 2nd step
    topPowers = simPowerTest_powerRef_2(length(simPowerTest_freqs)-i+1,:,:);
    topInacts = simPowerTest_inact_2(length(simPowerTest_freqs)-i+1,:,:);
    
    bottomPowers = simPowerTest_powerRef_2(i,:,:);
    bottomInacts = simPowerTest_inact_2(i,:,:);
    
    meanInact = zeros(1,length(simPowerTest_powers));
    
    for j = 1:length(simPowerTest_powers)
        meanInact(j) = mean([permute(topInacts(1,j,:), [2 3 1]) permute(bottomInacts(1,j,:), [2 3 1])]);
        
        plot(log10(permute(topPowers(1,j,:), [3 2 1])), permute(topInacts(1,j,:), [3 2 1])-meanInact(j), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'r');
        
        plot(log10(permute(bottomPowers(1,j,:), [3 2 1])), permute(bottomInacts(1,j,:), [3 2 1])-meanInact(j), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'b');
    end
    
    plot(log10(permute(topPowers(1,:,1), [3 2 1])), mean(permute(topInacts(1,:,:), [3 2 1]))-meanInact, '-.', 'linewidth', 2, 'markersize', 4, 'color', 'r')
    
    plot(log10(permute(bottomPowers(1,:,1), [3 2 1])), mean(permute(bottomInacts(1,:,:), [3 2 1]))-meanInact, '-.', 'linewidth', 2, 'markersize', 4, 'color', 'b')
    
    line([1 1]*log10(powerThresholdInterp_powerTest(i)), [-10 10], 'color', 'k')
end