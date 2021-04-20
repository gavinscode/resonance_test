compare_influenza_data

%% - look at influenza distributions
influenzaSize_mean = 100*10^-9;
influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

%%% Should adjust this based on freq found

influenzaVl = 1486; % Just matched in Saviot tool
influenzaVt = 743; %given x2 ratio

influenzaSize_samples = randn(500,1)*influenzaSize_std + influenzaSize_mean;

skewness(influenzaSize_samples)
kurtosis(influenzaSize_samples)

[influenzaSize_dist, influenzaSize_x] = hist(influenzaSize_samples,20);

influenzSize_resonances = zeros(length(influenzaSize_x),1);

tempResonance = [];

for i = 1:length(influenzaSize_x)
    influenzSize_resonances(i) = calcualtesphereresonance(influenzaSize_x(i)/2, ...
                'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
            
    tempResonance = [tempResonance, influenzSize_resonances(i)*ones(1,length(influenzaSize_dist(i)))];        
end

% Skewness increases and kurtosis decreases with freq transform..
skewness(tempResonance)
kurtosis(tempResonance)

figure;
subplot(1,2,1);
plot(influenzaSize_x*10^9, influenzaSize_dist)

subplot(1,2,2);
plot(influenzSize_resonances/10^9, influenzaSize_dist)


%% second plan - fitting curves for power and point for time

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
% Better to set frequency test at 75, 50, 25 etc.
simPowerTest_freqSpacing = 2;
% not sure of best spacing to use hear...
simPowerTest_powers = round(10.^([1:0.1:2, 2.5]))
simPowerTest_time = simFreqTest_time;

powerCols = cool(length(simPowerTest_powers));

% Interpolated to an s-curve
%%% would be simpler to just do with a sigmoid to keep things analytical
% Start with flat cut off
% Skew should be placed across range of size distribution, not inactivaiton
powerThresholdFreqs = [3 6 8.5 8.5+2.5 8.5+5.5]; % update if spread changed
powerThreshold = 60*[1 1 1 1 1]; 60*[1.5 1.3 0.9 0.7 0.4];

useWeibullPower = 1;

% Manual switch between ... 
% For full valued expression
powerWeibullAlpha = 0.008; % Scale
powerWeibullBeta = 0.9; % Shape

%%% Linear or Rayleigh CDF would also work well on suns data

% For log 10 expression
powerWeibullAlpha = 1.5; % Scale %2.5
powerWeibullBeta = 1.75; % Shape %1.75

% Linear expression - note this is log10 of power
powerLinearA = 0.709;
powerLinearB = -1.08;

powerRootExp = 0.4;
powerRootMax = 810; %%% Set so that curve hits max value around influenza

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
powerThresholdInterp_freqTest = interp1(powerThresholdFreqs, powerThreshold, simFreqTest_freqs, 'linear',0);
powerThresholdInterp_freqTest(powerThresholdInterp_freqTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
    simFreqTest_freqs(powerThresholdInterp_freqTest == 0), 'nearest','extrap');

for i = 1:length(simFreqTest_freqs)

    curveVal = curveMax*exp(-(simFreqTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    % Will always be above on this phase, but keep for consistancy
    if simFreqTest_power > powerThresholdInterp_freqTest(i)
%         if useWeibullPower
%             powerVal = (1-exp(-powerWeibullAlpha * ...
%                 ((log10(simFreqTest_power)-log10(powerThresholdInterp_freqTest(i))))^powerWeibullBeta));
%         else
%             powerVal = ((simFreqTest_power-powerThresholdInterp_freqTest(i))/powerRootMax)^powerRootExp; % 
%         end
        
        powerVal = powerLinearA*(log10(simFreqTest_power-powerThresholdInterp_freqTest(i)))+powerLinearB;
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
freq_curveConfidence = predint(freqCurve, freqRangeFine, 0.95, 'Functional');
freq_curveObserved = predint(freqCurve, freqRangeFine, 0.95, 'obs');

% Get inactivation limits
tempInds = find(freq_curveConfidence(:,1) > simPowerTest_inactLimit);
freqLimitsInact = freq_curveConfidence(tempInds([1 end]),1);
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
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'rd', 'linewidth', 2, 'markersize', 4)
xlim([1 20]); ylim([-10 110])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);
% plot(freqRangeFine, freq_curveObserved, 'r:', 'linewidth', 2)

% plot(simFreqTest_FreqRef(:,1), simFreqTest_inactRef(:,1), 'm-')

% plot(tempFreqs(1), freqLimitsInact(1), 'cx', 'markersize', 8, 'linewidth', 2);
% plot(tempFreqs(2), freqLimitsInact(2), 'cx', 'markersize', 8, 'linewidth', 2);

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

%%% Test deconvolution - difference ends up basically <1% - too small...
subplot(1,3,3); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2)/max(dataInactFreq2016(:,2))*freqCurve.a1, '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
% xlim([1 20]); ylim([-10 110])

freqCurve_scaled = freqCurve;

freq_curveInact_scaled = freqCurve_scaled(freqRangeFine);

% Interpolate influenza to same spacing
influenzaSize_res_fine = freqRangeFine; 
influenzaSize_res_fine(influenzaSize_res_fine < (min(influenzSize_resonances)) |...
    influenzaSize_res_fine > (max(influenzSize_resonances))) = [];
influenzaDist_fine = interp1(influenzSize_resonances, influenzaSize_dist, influenzaSize_res_fine, 'linear'); 

% Interpolate absorbtion to same spacing
absorbtion_res_fine = freqRangeFine;
absorbtion_res_fine(absorbtion_res_fine < min(dataAbs2016(:,1)) |...
    absorbtion_res_fine > max(dataAbs2016(:,1))) = [];
absorbtion_fine = interp1(dataAbs2016(:,1), dataAbs2016(:,2), absorbtion_res_fine, 'linear'); 

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);

% inactFn = conv(freq_curveInact, influenzaDist_fine,'same')/sum(influenzaDist_fine);
% inactFn = deconv(freq_curveInact_scaled, influenzaDist_fine);
% plot(freqRangeFine, inactFn, 'm-')

inactFn = conv(absorbtion_fine, influenzaDist_fine,'full')/sum(influenzaDist_fine);
plot(freqRangeFine(90:(90+length(inactFn)-1)), inactFn/max(inactFn)*freqCurve.a1, 'm-')

plot(influenzaSize_res_fine, influenzaDist_fine, 'b-')

plot(absorbtion_res_fine, absorbtion_fine/max(absorbtion_fine)*freqCurve.a1, 'k-')

% difference shows conv has little effect.
% figure;
% plot(freq_curveInact - inactFn)

%% Phase 1.2 - power
simPowerTest_inact = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_inactRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_freqRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);
simPowerTest_powerRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps);

% Get matched roots across freq
powerThresholdInterp_powerTest = interp1(powerThresholdFreqs, powerThreshold, simPowerTest_freqs, 'linear',0);
powerThresholdInterp_powerTest(powerThresholdInterp_powerTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
    simPowerTest_freqs(powerThresholdInterp_powerTest == 0), 'nearest','extrap');

for i = 1:length(simPowerTest_freqs)

    curveVal = curveMax*exp(-(simPowerTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simPowerTest_powers)

        if simPowerTest_powers(j) > powerThresholdInterp_powerTest(i)
%             if useWeibullPower
%                 powerVal = (1-exp(-powerWeibullAlpha * ...
%                     ((log10(simPowerTest_powers(j))-log10(powerThresholdInterp_powerTest(i))))^powerWeibullBeta));
%             else
%                 powerVal = ((simPowerTest_powers(j)-powerThresholdInterp_powerTest(i))/powerRootMax)^powerRootExp; % 
%             end

            powerVal = powerLinearA*(log10(simPowerTest_powers(j)-powerThresholdInterp_powerTest(i)))+powerLinearB;
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

%%% Maybe just use this where freq is between power refs
freqsToUse = find(simFreqTest_FreqRef >= simPowerTest_freqs(1) &  simFreqTest_FreqRef <= simPowerTest_freqs(end));
[testFreq, testPower] = meshgrid(freqRange, powerRange);

fit2D = scatteredInterpolant([simPowerTest_freqRef(:)' simFreqTest_FreqRef(freqsToUse)']', log10([simPowerTest_powerRef(:)', simFreqTest_PowerRef(freqsToUse)']'),...
    [simPowerTest_inact(:)', simFreqTest_inact(freqsToUse)']','linear','none');

interpData = fit2D(testFreq, log10(testPower));
interpData(isnan(interpData(:))) = 0;
% 
% fit2D = fit([[simPowerTest_freqRef(:)' simFreqTest_FreqRef(freqsToUse)']', log10([simPowerTest_powerRef(:)', simFreqTest_PowerRef(freqsToUse)']')],...
%     [simPowerTest_inact(:)', simFreqTest_inact(freqsToUse)']','lowess');
% 
% interpData = fit2D(testFreq, log10(testPower));

% try gaussian kernel
% interpData = ksrmv([simPowerTest_freqRef_2, log10(simPowerTest_powerRef_2)], simPowerTest_inact_2); %, [], [testFreq(:), testPower(:)]);
% interpData.h
% interpData = ksrmv([simPowerTest_freqRef_2, log10(simPowerTest_powerRef_2)], simPowerTest_inact_2,...
%     interpData.h/4, [testFreq(:), log10(testPower(:))]);
% 
% interpData = reshape(interpData.f, size(testFreq));

%%% Get sigma from original... 
%%% Add middle curve in.
% gkrModel = fitrgp([[simPowerTest_freqRef(:)' simFreqTest_FreqRef(freqsToUse)']', log10([simPowerTest_powerRef(:)', simFreqTest_PowerRef(freqsToUse)']')], ...
%     [simPowerTest_inact(:)', simFreqTest_inact(freqsToUse)']')
% 
% [interpData,ysd,yint] = predict(gkrModel, [testFreq(:), log10(testPower(:))]);


figure;
subplot(1,3,2)
plotData = reshape(interpData, size(testFreq));
imshow(plotData/100); hold on
pause(0.1)

lowPoints = find(plotData < 1);
[lowX, lowY] = ind2sub(size(plotData), lowPoints);

plot(lowY, lowX, '.r')

freqRef = zeros(length(simPowerTest_freqs),1);
powerRef = zeros(length(simPowerTest_freqs),1);

for i = 1:length(simPowerTest_freqs)
   [~, freqRef(i)] = min(abs(freqRange - simPowerTest_freqs(i)));
   [~, powerRef(i)] = min(abs(log10(powerRange) - log10(powerThresholdInterp_powerTest(i))));
end

plot(freqRef, powerRef, 'm-o')

% plotData = reshape(yint(:,1), size(testFreq));
subplot(1,3,1)
imshow(plotData/100); hold on
pause(0.1)

lowPoints = find(plotData < 2);
[lowX, lowY] = ind2sub(size(plotData), lowPoints);

plot(lowY, lowX, '.r')
plot(freqRef, powerRef, 'm-o')

% plotData = reshape(yint(:,2), size(testFreq));
subplot(1,3,3)
imshow(plotData/100); hold on
pause(0.1)

lowPoints = find(plotData < 5);
[lowX, lowY] = ind2sub(size(plotData), lowPoints);

plot(lowY, lowX, '.r')
plot(freqRef, powerRef, 'm-o')

figure; 

opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 freqCurve.c1/2], 'Upper', [freqCurve.a1 freqCurve.b1 freqCurve.c1]);

for i = fliplr(1:length(simPowerTest_powers))
    
    subplot(3,7,1); hold on
    
    tempFreqs = simPowerTest_freqRef(:,i,:);
    tempInacts = simPowerTest_inact(:,i,:);
    
    plot(permute(tempFreqs, [1 3 2]), permute(tempInacts, [1 3 2]), 'o', 'linewidth', 2, 'markersize', 4, 'color', powerCols(i,:))

    if i ~= length(simPowerTest_powers)
        % Take limits from precedding
        
        opts = fitoptions('gauss1', 'Lower', [0 freqCurve.b1 powerCurve.c1/2], 'Upper', [powerCurve.a1 freqCurve.b1 powerCurve.c1]);
    end
        
    powerCurve = fit(tempFreqs(:), tempInacts(:), 'gauss1', opts)

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
    
    % get freq bounds
    subplot(3,7,2); hold on
    
    tempInds = find(power_curveInact > simPowerTest_inactLimit);
    if ~isempty(tempInds)
        tempFreqs = [freqRangeFine(tempInds(1)) freqRangeFine(tempInds(end))];
    end
    
    %%% Sometimes these end up being very narrow - not sure if that is a problem or emergent feature.
    tempFreqsInner = []; 
    tempFreqsOuter = [];
    
    if ~all( all( isnan( confint(powerCurve))))
        % Get outer
        tempInds = find(power_curveConfidence(:,1) > simPowerTest_inactLimit);
        if ~isempty(tempInds)
           tempFreqsInner = [freqRangeFine(tempInds(1)) freqRangeFine(tempInds(end))];
        end
        
        % Get inner
        tempInds = find(power_curveConfidence(:,2) > simPowerTest_inactLimit);
        if ~isempty(tempInds)
           tempFreqsOuter = [freqRangeFine(tempInds(1)) freqRangeFine(tempInds(end))];
        end
    end
    
    if ~isempty(tempFreqsOuter) & ~isempty(tempFreqsInner)
        line([tempFreqsOuter(1) tempFreqsInner(1)], [1 1]*simPowerTest_powers(i), 'color', powerCols(i,:), 'linewidth', 1)
        line([tempFreqsInner(2) tempFreqsOuter(2)], [1 1]*simPowerTest_powers(i), 'color', powerCols(i,:), 'linewidth', 1)
        
        plot(tempFreqs(1), simPowerTest_powers(i), 'o', 'markersize', 8, 'linewidth', 1, 'color', powerCols(i,:));
        plot(tempFreqs(2), simPowerTest_powers(i), 'o', 'markersize', 8, 'linewidth', 1, 'color', powerCols(i,:));
    end
    
    % Test convultion
    subplot(3,7,3); hold on
    inactFn = conv(power_curveInact, influenzaDist_fine,'same')/sum(influenzaDist_fine);
%     plot(power_curveInact - inactFn, 'color', powerCols(i,:))
    plot(inactFn, 'color', powerCols(i,:))
%     inactFn = deconv(power_curveInact, influenzaDist_fine);
%     plot(inactFn, 'color', powerCols(i,:))
    
    subplot(3,7,4); hold on
    plot(simPowerTest_freqRef(:,i,1), simPowerTest_inactRef(:,i,1), 'color', powerCols(i,:));
    
    subplot(3,7,5); hold on
    plot(freqRangeFine, power_curveInact, 'color', powerCols(i,:));
end

subplot(3,7,1); hold on

plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2) - (100-freqCurve.a1), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
xlim([0 20]); ylim([0 100])

subplot(3,7,2); hold on
plot(simPowerTest_freqs, powerThresholdInterp_powerTest, 'k')

plot(influenzSize_resonances/10^9, influenzaSize_dist+200, '-b', 'linewidth', 2)

for i = 1:length(simPowerTest_freqs)
    
    tempPowers = simPowerTest_powerRef(i,:,:);
    tempInacts = simPowerTest_inact(i,:,:);

    
    % Plot log
    subplot(3,7,7+i); hold on
    plot(log10(permute(tempPowers, [2 3 1])), permute(tempInacts, [2 3 1]), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'r')

    if i == centrePowerInd
        plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    end
    
    % plot normal
    subplot(3,7,14+i); hold on
    plot((permute(tempPowers, [2 3 1])), permute(tempInacts, [2 3 1]), 'o', 'linewidth', 2, 'markersize', 4, 'color', 'r')

    if i == centrePowerInd
        plot((dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    end
    
    % Do fit
    offSetValue = topCurve(simPowerTest_freqs(i));
    
    weibullFit=fittype('(x>c)*(1-exp(-a*((x-c))^b))',...
             'coefficients',{'a','b','c'},'independent','x');
    
    opts = fitoptions(weibullFit);
    opts.Lower = [0.5 0.5 0];
    opts.Upper = [5 4 log10(100)];
    opts.StartPoint = [1 1 log10(50)];
    
    fittedCurve = fit(log10(tempPowers(:)), tempInacts(:)/offSetValue, weibullFit, opts);
    
    % offSetValue = 1;
%     linearFit = fittype('(a*(x)+b)',... %  * ((a*(x-c)+b)<100) + ((a*(x-c)+b)>100)*100
%              'coefficients',{'a','b'},'independent','x');
%     
%     opts = fitoptions(linearFit);
% %     opts.Lower = [0 -20 0];
% %     opts.Upper = [10 20 log10(100)];
%     opts.StartPoint = [1 0];
%     
%     fittedCurve = fit(log10(tempPowers(:)), tempInacts(:), linearFit, opts);

    power_curveInact = fittedCurve(log10(powerRangeFine))*offSetValue;
    power_coefBound = confint(fittedCurve);
    
    % Plot log
    subplot(3,7,7+i); hold on
    plot(log10(powerRangeFine), power_curveInact, 'color', 'r')
    
    if ~all( all( isnan( power_coefBound)))
        power_curveConfidence = predint(fittedCurve, log10(powerRangeFine), 0.95, 'Functional')*offSetValue;
        power_curveObserved = predint(fittedCurve, log10(powerRangeFine), 0.95, 'obs')*offSetValue;
        
        plot(log10(powerRangeFine), power_curveConfidence, '--', 'color', 'r')
    else
        power_curveConfidence = 0;
        power_curveObserved = 0;
    end
    
    xlim([1 3]); ylim([0 100])
    
    if power_coefBound(2,3) < 3
        title(sprintf('%.1f %.1f:%.1f:%.1f', powerThresholdInterp_powerTest(i), 10^power_coefBound(1,3), ...
            10^fittedCurve.c, 10^power_coefBound(2,3)));
    else
        title(sprintf('%.1f', powerThresholdInterp_powerTest(i)));
    end
    % Plot ref curve
    powersToUse = powerRangeFine(powerRangeFine > powerThresholdInterp_powerTest(i));
    
    curveVal = curveMax*exp(-(simPowerTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
%     powerVal = (1-exp(-powerWeibullAlpha .* ...
%                     ((log10(powersToUse)-log10(powerThresholdInterp_powerTest(i)))).^powerWeibullBeta));
    powerVal = powerLinearA*(log10(powersToUse-powerThresholdInterp_powerTest(i)))+powerLinearB;           
    totalVal = curveVal*powerVal;            
    totalVal(totalVal > 100) = 100; 
    totalVal(totalVal < 0) = 0;
    
    plot(log10(powersToUse), totalVal, 'k')
    
    % Plot normal
    subplot(3,7,14+i); hold on
    plot((powerRangeFine), power_curveInact, 'color', 'r')
    
    if ~all( all( isnan( power_coefBound)))
        plot((powerRangeFine), power_curveConfidence, '--', 'color', 'r')
    end
    
    xlim([0 1000]); ylim([0 100])

    plot((powersToUse), curveVal*powerVal, 'k')
    
    
    % get power bounds
    subplot(3,7,2); hold on
    tempInds = find(power_curveInact > simPowerTest_inactLimit);
    if ~isempty(tempInds)
        tempPower = powerRangeFine(tempInds(1));
    end
    
    tempPowerInner = [];
    tempPowerOuter = [];
    
    if ~all( all( isnan( power_coefBound)))
        % Get outer
        tempInds = find(power_curveConfidence(:,1) > simPowerTest_inactLimit);
        if ~isempty(tempInds)
           tempPowerInner = powerRangeFine(tempInds(1));
        end
        
        % Get inner
        tempInds = find(power_curveConfidence(:,2) > simPowerTest_inactLimit);
        if ~isempty(tempInds)
           tempPowerOuter = powerRangeFine(tempInds(1));
        end
    end
    
    if ~isempty(tempPowerOuter) & ~isempty(tempPowerInner)
        line([1 1]*simPowerTest_freqs(i), [tempPowerInner tempPowerOuter], 'linewidth', 1)
        
        plot(simPowerTest_freqs(i), tempPower, 'o', 'markersize', 8, 'linewidth', 1);
        
        subplot(3,7,7+i)
%         title(sprintf('%.1f %.1f', powerThresholdInterp_powerTest(i), tempPower));
    end
end