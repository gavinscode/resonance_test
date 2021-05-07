compare_influenza_data
close all
%% - generate parameters for size - would be from SEM

influenzaSize_mean = 100*10^-9;
influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

initialCountNum = 5000;

influenzaSize_samples = randn(initialCountNum,1)*influenzaSize_std + influenzaSize_mean;

% How accurately can size be determined from TEM? surely not < 1 nm... 

distStep = 1; % Also percentage as across 100
stepPercent = distStep;
%%% May be some optimum way to set step to get best SNR later.

[influenzaSize_dist, influenzaSize] = hist(influenzaSize_samples,(50:distStep:150)/10^9);

[~, maxSizeInd] = max(influenzaSize_dist);

%%% This should actually be done after we know peak freq and largest size
influenzaVl = 1486; % Just matched in Saviot tool
influenzaVt = 743; % given x2 ratio

tempResSmall = calcualtesphereresonance(influenzaSize(1)/2, ...
                'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
            
tempResLarge = calcualtesphereresonance(influenzaSize(end)/2, ...
    'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;

% Note resonance of small is higher than large
resSlope = (tempResSmall-tempResLarge)/(1/(influenzaSize(1)/2)-1/(influenzaSize(end)/2));

influenzSize_resonances = 1./(influenzaSize/2)*resSlope;

% influenzSize_resonances = zeros(length(influenzaSize),1);

% for i = 1:length(influenzaSize)
%     influenzSize_resonances(i) = calcualtesphereresonance(influenzaSize(i)/2, ...
%                 'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
% end

measuredSize_resonances = 1./(dataDiameter2017(:,1)/2/10^9)*resSlope;

figure;
subplot(1,2,1); hold on
plot(influenzaSize, influenzaSize_dist);

plot(dataDiameter2017(:,1)/10^9, dataDiameter2017(:,2));

subplot(1,2,2); hold on
plot(influenzSize_resonances, influenzaSize_dist);

plot(measuredSize_resonances, dataDiameter2017(:,2));

%% forth plan - distribution based

absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% phenomenlogical parameters
% For freq
curveMax = 100;
curveCenter = 8.5; 
curveSpread = 2.3;

% For power
powerThreshold = 0; %45

% Linear expression - note this is on log10 of power
powerLinearA = 0.709;
powerLinearB = -1.08;
zeroInactPoint = 10^(-powerLinearB/powerLinearA)

nReps = 5;
testCountNum = initialCountNum;

% Phase 1.1 - also want to fit function to get peak 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 630;
simFreqTest_time = 1000;

% Phase 1.2
% freqs will be on center
% powerLogRange = 1:0.25:2.25;
powerLogRange = 1.5:0.05:1.75;
simPowerTest_powers = round(10.^(powerLogRange))
simPowerTest_time = simFreqTest_time;

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

        powerVal = powerLinearA*(log10(simFreqTest_power))+powerLinearB;

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

simPowerTest_freqs = freqCurve.b1;

% Look at distributions - take for first result
figure; hold on

%%% Make a function for this, will use later
for i = 1:length(simFreqTest_freqs)
    
    for j = 1:nReps
        tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

        tempSize_samples = sort(tempSize_samples, 'descend');

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

        samplesIntact = ones(testCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(testCountNum*simFreqTest_inactRef(i,1)/100);
        
        %%% Set flag 
%         numToInact = round(testCountNum*simFreqTest_inact(i,j)/100);

        trueInact = round(testCountNum*simFreqTest_inactRef(i,1)/100);
        
        if numToInact > 1 & trueInact > 1
            % Find nearest in sample
            [~, minInd] = min(abs(tempSize_freqs-simFreqTest_freqs(i)));

            numInactivated = 1;

            samplesIntact(minInd) = 0;

            while numInactivated < numToInact
                % Find nearest
                if minInd + 1 > testCountNum
                    minInd = minInd - 1;
                elseif minInd - 1 < 1
                    minInd = minInd + 1;    
                else
                    indRef = find(samplesIntact);

                    [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-simFreqTest_freqs(i)));

                    minInd = indRef(tempInd);
                end

                samplesIntact(minInd) = 0;

                numInactivated = numInactivated + 1;
            end
        end

        % Get hist of both
        inactiveDist = hist(tempSize_samples(samplesIntact == 0),influenzaSize);

        activeDistDist = hist(tempSize_samples(samplesIntact == 1),influenzaSize);

        subplot(2,10,i); hold on
        plot(influenzaSize, activeDistDist+inactiveDist, 'r')

        plot(influenzaSize, activeDistDist, 'k')
    end
end

figure;
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'rd', 'linewidth', 2, 'markersize', 4)
xlim([1 20]); ylim([-10 110])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);
% plot(freqRangeFine, freq_curveObserved, 'r:', 'linewidth', 2)

% plot(simFreqTest_FreqRef(:,1), simFreqTest_inactRef(:,1), 'm-')

% Show resonance distribution
plot(influenzSize_resonances, influenzaSize_dist, '-b', 'linewidth', 2)

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

% For linear intersect to zero
% simPowerTest_powers(1:5) = 33:0.25:34;

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
            powerVal = powerLinearA*(log10(simPowerTest_powers(j)))+powerLinearB;    
            
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

% Plot power curve
figure; 
subplot(3,length(simPowerTest_powers),1); hold on

for j = 1:length(simPowerTest_powers)
    plot(log10(permute(simPowerTest_powerRef(1,j,:), [3 2 1])), permute(simPowerTest_inact(1,j,:), [3 2 1]), 'ro', 'linewidth', 2, 'markersize', 4)
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])


inactRatioBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps, length(influenzaSize));

inactFractionBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps, length(influenzaSize));

activeBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nReps, length(influenzaSize));

%%% Main to do here
%%%     1. scale active virus particle dist sides to match original dist
%%%     2. Pick drop in dist/peak innact and track as it gets lower - otherwise centre point needs to be regularly adjusted
            % May move by dist

%%% Make a function for this, will use later
for i = 1:length(simPowerTest_powers)
    
    subplot(3,length(simPowerTest_powers),length(simPowerTest_powers)+i); hold on
    plot(influenzaSize_dist/initialCountNum*100, '-b', 'linewidth', 2)
        
    for j = 1:nReps
        tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

        tempSize_samples = sort(tempSize_samples, 'descend');

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

        samplesIntact = ones(testCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(testCountNum*simPowerTest_inactRef(1,i,1)/100);
        
        %%% Set flag 
%         numToInact = round(testCountNum*simPowerTest_inact(1,i,j)/100);
        trueInact = round(testCountNum*simPowerTest_inactRef(1,i,1)/100);
        
        if numToInact > 1 & trueInact > 1
            % Find nearest in sample
            [~, minInd] = min(abs(tempSize_freqs-simPowerTest_freqs(1)));

            numInactivated = 1;

            samplesIntact(minInd) = 0;

            while numInactivated < numToInact
                % Find nearest
                if minInd + 1 > testCountNum
                    minInd = minInd - 1;
                elseif minInd - 1 < 1
                    minInd = minInd + 1;    
                else
                    indRef = find(samplesIntact);

                    [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-simPowerTest_freqs(1)));

                    minInd = indRef(tempInd);
                end

                samplesIntact(minInd) = 0;

                numInactivated = numInactivated + 1;
            end
        end

        % Get hist of both
        inactiveDist = hist(tempSize_samples(samplesIntact == 0),influenzaSize);

        activeDistDist = hist(tempSize_samples(samplesIntact == 1),influenzaSize);

        % Strictly speaking, to match equivelent of usual Ct/C0 doesn't require distribution of inactivated
        %%% However, should know proper ratio, which may require fitting sides of active distribution
        zeroRatioInds = find(influenzaSize_dist/initialCountNum*100 < stepPercent);
        
        zeroFractionInds = find((activeDistDist+inactiveDist)/testCountNum*100 < stepPercent);
        
        inactRatioBySize(1,i,j,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;
        
        inactRatioBySize(1,i,j,zeroRatioInds) = 0;
        
        inactFractionBySize(1,i,j,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;
        
        inactFractionBySize(1,i,j,zeroFractionInds) = 0;
        
        activeBySize(1,i,j,:) = activeDistDist/testCountNum*100;
        
        % Counted inactivated
        subplot(3,length(simPowerTest_powers),length(simPowerTest_powers)+i); hold on

        plot((activeDistDist+inactiveDist)/testCountNum*100, 'r') %activeDistDist+

%         plot(influenzaSize, (inactiveDist)/testCountNum*100, 'r')
        
        plot(activeDistDist/testCountNum*100, 'k')
        
        %inactivation spectrum
        subplot(3,length(simPowerTest_powers),2*length(simPowerTest_powers)+i); hold on

        plot(permute(inactRatioBySize(1,i,j,:), [4 3 2 1]), 'k');
        ylim([-50 100])
    end
end

% Getting this spot on is hard - may need to do manually...
% Need to set this better
if distStep == 1
    %%% Somehow this shifts around depending on number of samples...
    % moves by dist
    maxFreqInd = 52; % for 1 step...
    %51 for 1000, 52 for 1000
elseif distStep == 0.5
    maxFreqInd = 101; % for 0.5 step...
    % seems more stable against shifting
else
    error('Need to find min diststep')
end

% Given count
subplot(3,length(simPowerTest_powers),2); hold on

for j = 1:length(simPowerTest_powers)
    plot(log10(permute(simPowerTest_powerRef(1,j,:), [3 2 1])), permute(inactFractionBySize(1,j,:, maxFreqInd), [4, 3 2 1]), 'ro', 'linewidth', 2, 'markersize', 4)
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])

% Given original
subplot(3,length(simPowerTest_powers),3); hold on

for j = 1:length(simPowerTest_powers)
    plot(log10(permute(simPowerTest_powerRef(1,j,:), [3 2 1])), permute(inactRatioBySize(1,j,:,maxFreqInd), [4 3 2 1]), 'ro', 'linewidth', 2, 'markersize', 4)
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])

