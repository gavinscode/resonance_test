compare_influenza_data

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
 
%% second plan - sequential power steps.

nReps = 3;
absStd= 10; % std on absolute, not relative values (lower SNR on low inactiviation)

% Phase 1.1 - also want to fit function to get peak 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 630;
simFreqTest_time = 1000;

curveMax = 100;
curveCenter = 8.5; 
curveSpread = 2.3;

% Phase 1.2
simPowerTest_inactLimit = 1; % from height of prediction bound, will be on min and max
% 1st freq points from inact distribution, 2nd from size dist. 
    % Negative is below peak
simPowerTest_freqPoints = {[1, 0.5, -0.5], [0.5 -0.5]};
% not sure of best spacing to use here...
simPowerTest_powers = round(10.^(1:0.25:2.25));
simPowerTest_time = simFreqTest_time;

powerCols = cool(length(simPowerTest_powers));

% Just use single power threshold
powerThreshold = 56;

useSharpStep = 1;
% need to adjust curve if not...

% Power decay
% Linear expression - note this is on log10 of power
powerLinearA = 0.709;
powerLinearB = -1.08;

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

% Get matched roots across freq
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
   end
end

% Get freq curve from size with match spacing
resonances_fine = interp1(influenzSize_resonances, influenzaSize_dist, freqRangeFine, 'linear', 0);
resonances_fine(isnan(resonances_fine)) = 0;

sizeCurveRatios = simPowerTest_freqPoints{2};
sizeCurveFreqs = zeros(length(sizeCurveRatios), 1);

for i = 1:length(sizeCurveRatios)
   if sizeCurveRatios(i) ~= 1
       freqInd = find(resonances_fine/max(resonances_fine) > abs(sizeCurveRatios(i)));
       
       if sizeCurveRatios(i) > 0
           sizeCurveFreqs(i) = freqRangeFine(freqInd(end));
       else
           sizeCurveFreqs(i) = freqRangeFine(freqInd(1));
       end
   else
       error('1 should not be in size points')
   end
end

simPowerTest_freqs = [inactCurveFreqs' sizeCurveFreqs'];

simPowerTest_freqs = sort(simPowerTest_freqs);

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

opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 freqCurve.c1/2], 'Upper', [freqCurve.a1 freqCurve.b1 freqCurve.c1]);

for i = fliplr(1:length(simPowerTest_powers))
    
    subplot(1,3,1); hold on
    
    tempFreqs = simPowerTest_freqRef(:,i,:);
    tempInacts = simPowerTest_inact(:,i,:);
    
    plot(permute(tempFreqs, [1 3 2]), permute(tempInacts, [1 3 2]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols(i,:))

    if i ~= length(simPowerTest_powers)
        % Take limits from precedding
        
        opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 powerCurve.c1/2], 'Upper', [powerCurve.a1 freqCurve.b1 powerCurve.c1]);
    end
        
    powerCurve = fit(tempFreqs(:), tempInacts(:), 'gauss1', opts);
    power_coefBound = confint(powerCurve);

    if i == length(simPowerTest_powers)
       topCurve =  powerCurve;
    end
    
    power_curveInact = powerCurve(freqRangeFine);
    
    plot(freqRangeFine, power_curveInact, 'color', powerCols(i,:))
        
     if ~all( all( isnan( confint(powerCurve))))
        power_curveConfidence = predint(powerCurve, freqRangeFine, 0.95, 'Functional');
        power_curveObserved = predint(powerCurve, freqRangeFine, 0.95, 'obs');
        
%         plot(freqRangeFine, power_curveConfidence, '--', 'color', powerCols(i,:))
     else
         power_curveConfidence = [];
         power_curveObserved = [];    
     end
    
    % plot heigh
    subplot(1,3,2); hold on
    plot(log10(simPowerTest_powers(i)), powerCurve.a1, 'o', 'color', powerCols(i,:));
    line([1 1]*log10(simPowerTest_powers(i)), [power_coefBound(1,1), power_coefBound(2,1)], 'color', powerCols(i,:))
    
%     subplot(3,7,4); hold on
%     plot(simPowerTest_freqRef(:,i,1), simPowerTest_inactRef(:,i,1), 'color', powerCols(i,:));
    
end

subplot(1,3,1); hold on

% plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
xlim([0 20]); ylim([-10 110])

subplot(1,3,2); hold on
line([1 1]*log10(powerThreshold), [0 100], 'color', 'k')

plot(log10(simFreqTest_power(i)), freqCurve.a1, 'o', 'color', 'r');
line([1 1]*log10(simFreqTest_power(i)), [freq_coefBound(1,1), freq_coefBound(2,1)], 'color', 'r')

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])

