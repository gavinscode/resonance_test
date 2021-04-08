compare_influenza_data

%%% second plan - fitting curves for power and point for time

nReps = 3;
absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% Phase 1.1 - also want to fit function to get peak 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 750;
simFreqTest_time = 1000;

curveMax = 100;
curveCenter = 8.5; 
curveSpread = 2.3;

% Phase 1.2
simPowerTest_inactLimit = 1; % from height of prediction bound, will be on min and max
simPowerTest_freqSpacing = 2;
simPowerTest_powers = [10 25 50 75 100 200 400];
simPowerTest_time = simFreqTest_time;

powerCols = cool(length(simPowerTest_powers));

% Interpolated to an s-curve
%%% would be simpler to just do with a sigmoid to keep things analytical
powerRootFreqs = [3 6 8.5 8.5+2.5 8.5+5.5]; % update if spread changed
powerRootMinimum = 55*[1.5 1.3 0.9 0.7 0.4];

powerRootExp = 0.4;
powerRootMax = 810; %%% Set so that curve hits max value around influenza

% Make simulated plots
% Set up image
freqRange = min(simFreqTest_freqs):0.2:max(simFreqTest_freqs);
freqRangeFine = 0:0.05:max(simFreqTest_freqs);
powerRange = 0:10:1000;

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
simFreqTest_FreqRef = zeros(length(simFreqTest_freqs), nReps);

% Get matched roots across freq
powerRootMin_matched_freqTest = interp1(powerRootFreqs, powerRootMinimum, simFreqTest_freqs, 'linear',0);
powerRootMin_matched_freqTest(powerRootMin_matched_freqTest == 0) = interp1(powerRootFreqs, powerRootMinimum, ...
    simFreqTest_freqs(powerRootMin_matched_freqTest == 0), 'nearest','extrap');

for i = 1:length(simFreqTest_freqs)

    curveVal = curveMax*exp(-(simFreqTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    % Will always be above on this phase, but keep for consistancy
    if simFreqTest_power > powerRootMin_matched_freqTest(i)
        simFreqTest_inact(i,:) = ((simFreqTest_power-powerRootMin_matched_freqTest(i))/powerRootMax)^powerRootExp *curveVal + ...
            absStd*randn(nReps,1); % 
    else
        simFreqTest_inact(i,:) = absStd*randn(nReps,1); %
    end

    simFreqTest_FreqRef(i,:) = simFreqTest_freqs(i);
end

% fit curve
freqCurve = fit(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'gauss1');

freq_curveInact = freqCurve(freqRangeFine);
freq_curveConfidence = predint(freqCurve, freqRangeFine, 0.95, 'Functional');
freq_curveObserved = predint(freqCurve, freqRangeFine, 0.95, 'obs');

% Get inactivation limits
tempInds = find(freq_curveConfidence(:,1) > simPowerTest_inactLimit);
freqLimitsPowers = freq_curveConfidence(tempInds([1 end]),1);
tempFreqs = [freqRangeFine(tempInds(1)) freqRangeFine(tempInds(end))];

simPowerTest_FreqLimits = [floor(tempFreqs(1)) ceil(tempFreqs(2))];
if simPowerTest_FreqLimits(1) < 1
    simPowerTest_FreqLimits(1)  = 1;
end

centreFreq = round(freqCurve.b1*10)/10;

[~, centreFreqInd_nearest] = min(abs(simFreqTest_freqs - centreFreq));

% Now place freqs
simPowerTest_freqs = centreFreq;

% Add lower
while (simPowerTest_freqs(1) - simPowerTest_freqSpacing) > simPowerTest_FreqLimits(1)
    simPowerTest_freqs = [(simPowerTest_freqs(1) - simPowerTest_freqSpacing) simPowerTest_freqs];
end

% Add upper
while (simPowerTest_freqs(end) + simPowerTest_freqSpacing) < simPowerTest_FreqLimits(2)
    simPowerTest_freqs = [simPowerTest_freqs (simPowerTest_freqs(end) + simPowerTest_freqSpacing)];
end

centrePowerInd = find(simPowerTest_freqs == centreFreq);

figure;
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'rd', 'linewidth', 2, 'markersize', 4)
xlim([1 20]); ylim([0 100])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);
% plot(freqRangeFine, freq_curveObserved, 'r:', 'linewidth', 2)

plot(tempFreqs(1), freqLimitsPowers(1), 'cx', 'markersize', 8, 'linewidth', 2);
plot(tempFreqs(2), freqLimitsPowers(2), 'cx', 'markersize', 8, 'linewidth', 2);

% show image
%Place freq inactivation in map
%%% note, could do with fitted curve as well...
interpInactFreq = interp1(simFreqTest_freqs, mean(simFreqTest_inact,2),...
    freqRange, 'linear', 0);

interpInactFreq(interpInactFreq < 2) = 2;

[~, powerInd] = min(abs(powerRange - simFreqTest_power));
inactImage(:, powerInd) = interpInactFreq/100;

freqRef = zeros(length(simFreqTest_freqs),1);
for i = 1:length(simFreqTest_freqs)
   [~, freqRef(i)] = min(abs(freqRange - simFreqTest_freqs(i)));
end

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
simPowerTest_freqRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_powerRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);

% Get matched roots across freq
powerRootMin_matched_powerTest = interp1(powerRootFreqs, powerRootMinimum, simPowerTest_freqs, 'linear',0);
powerRootMin_matched_powerTest(powerRootMin_matched_powerTest == 0) = interp1(powerRootFreqs, powerRootMinimum, ...
    simPowerTest_freqs(powerRootMin_matched_powerTest == 0), 'nearest','extrap');

for i = 1:length(simPowerTest_freqs)

    curveVal = curveMax*exp(-(simPowerTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simPowerTest_powers)

        if simPowerTest_powers(j) > powerRootMin_matched_powerTest(i)
            simPowerTest_inact(i,j,:) = ((simPowerTest_powers(j)-powerRootMin_matched_powerTest(i))/powerRootMax)^powerRootExp *curveVal + ...
                absStd*randn(nReps,1); % 
        else
            simPowerTest_inact(i,j,:) = absStd*randn(nReps,1); %
        end
        
        simPowerTest_freqRef(i,j,:) = simPowerTest_freqs(i);
        simPowerTest_powerRef(i,j,:) = simPowerTest_powers(j);
    end
end

figure; 
subplot(1,2,1); hold on

opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 freqCurve.c1/2], 'Upper', [freqCurve.a1 freqCurve.b1 freqCurve.c1*1.5]);

for i = fliplr(1:length(simPowerTest_powers))
    tempFreqs = simPowerTest_freqRef(:,i,:);
    tempInacts = simPowerTest_inact(:,i,:);
    
    plot(permute(tempFreqs, [1 3 2]), permute(tempInacts, [1 3 2]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols(i,:))

    if i ~= length(simPowerTest_powers)
        % Take limits from precedding
        opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 powerCurve.c1/2], 'Upper', [powerCurve.a1 freqCurve.b1 powerCurve.c1]);
    end
        
    powerCurve = fit(tempFreqs(:), tempInacts(:), 'gauss1', opts);

    power_curveInact = powerCurve(freqRangeFine);
    
    plot(freqRangeFine, power_curveInact, 'color', powerCols(i,:))
        
    if powerCurve.a1 > 10^-3
        power_curveConfidence = predint(powerCurve, freqRangeFine, 0.95, 'Functional');
        power_curveObserved = predint(powerCurve, freqRangeFine, 0.95, 'obs');
        
        plot(freqRangeFine, power_curveConfidence, '--', 'color', powerCols(i,:))
    end
end

plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
xlim([0 20]); ylim([0 100])

subplot(1,2,2); hold on

tempPowers = simPowerTest_powerRef(centrePowerInd,:,:);
tempInacts = simPowerTest_inact(centrePowerInd,:,:);

plot(permute(tempPowers, [2 3 1]), permute(tempInacts, [2 3 1]), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'r')
plot(simFreqTest_power, simFreqTest_inact(centreFreqInd_nearest,:), 'rd', 'linewidth', 2, 'markersize', 4)

plot(dataInactPower2016(:,1), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])

xlim([0 1000]); ylim([0 100])