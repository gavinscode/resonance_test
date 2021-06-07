compare_influenza_data
close all
%% - generate parameters for size - would be from SEM

influenzaSize_mean = 100*10^-9;
influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

initialCountNum = 5000;

influenzaSize_samples = randn(initialCountNum,1)*influenzaSize_std + influenzaSize_mean;

% How accurately can size be determined from TEM? surely not < 1 nm... 

distStep = 2;

distPercent = 0.5; % Note that nothing below stepPercent will get included in size maps
stepPercent = 2;

%%% May be some optimum way to set step to get best SNR later.

[influenzaSize_dist, influenzaSize] = hist(influenzaSize_samples,(70:distStep:130)/10^9);

% Trim distribution and samples
tempInd = find(influenzaSize_dist > sum(influenzaSize_dist)*distPercent/100);

%%% Note that bin centers are specified
influenzaSize_samples(influenzaSize_samples < mean(influenzaSize(tempInd(1)-1:tempInd(1))) | ...
    influenzaSize_samples > mean(influenzaSize(tempInd(end):tempInd(end)+1))) = [];

influenzaSize(influenzaSize_dist < sum(influenzaSize_dist)*distPercent/100) = [];
influenzaSize_dist(influenzaSize_dist < sum(influenzaSize_dist)*distPercent/100) = [];

initialCountNum = sum(influenzaSize_dist);

if length(influenzaSize_samples) ~= initialCountNum
   error('Problem triming') 
end

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
    % Later - need to add replication and calibration test on distribution

absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% phenomenlogical parameters
    % For freq
        curveMax = 100;
        curveCenter = 8.5; 
        curveSpread = 2.3;

    % For power - single threhsold
        powerThreshold = 45; %45
        % For power - linear interpolation between points then constant
%         powerThreshold = [35 55]; 
%         powerThresholdFreqs = [8 9];

        % Linear expression - note this is on log10 of power
        powerLinearA = 0.709;
        powerLinearB = -1.08;
        zeroInactPoint = 10^(-powerLinearB/powerLinearA)

        % Weibull expression
        powerWeibullThreshold = log10(zeroInactPoint);
        powerWeibullAlpha = 1.5; % Scale
        powerWeibullBeta = 1.75; % Shape

        % Which to use
        useWeibullPower = 0;

        % For time
        timeWeibullThreshold = log10(0.2);
        timeWeibullAlpha = 1.5; % Scale
        timeWeibullBeta = 1.75; % Shape
        
        timeWeibullPowerCenter = 2; % Center value power is scaled aorund
        timeWeibullPowerScale = 6;  % Magnitude of scale
        
% Test parameters
nReps = 3; % For plaque assay
nrepsCount = 1; % for EM - will be match to one plaque assay, but independent in simulation

countVec = 1:nrepsCount;

testCountNum = initialCountNum;

% Phase 1.1 - also want to fit function to get peak 
% 2 GHz seperation seems ok to fit curve to inact, but then have less info on individual inact 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = round(10^2.75);
simFreqTest_time = 1000;

freq_freqCols = jet(length(simFreqTest_freqs)); 

% Phase 1.2 - identify power response (minimum effective) at center freq
% powerLogRange = [1:0.25:2.25]; % coarse search
%%% Not neccersary to duplicate powers from previous step i.e. start at 1.55 and 1.7
% powerLogRange = 1.5:0.05:1.75; % For fine search. works on 45 or 33

powerLogRange = [1:0.25:2.75];
% simPowerTest_freqs = [6 8.5 11]; 
% simPowerTest_freqs = [7.5 8.5 9.5]; % Given 100 nm, 7 and 10 GHz are limits of inactivation (6 and 11 sometimes included) . halfway between these and center
simPowerTest_freqs = 8.5;

simPowerTest_powers = round(10.^(powerLogRange))
simPowerTest_time = simFreqTest_time;

%%% Need to do second step to get fine values

power_freqCols = winter(length(simPowerTest_freqs));
power_powerCols = cool(length(simPowerTest_powers));

% Phase 1.3 - scan across freq with lower powers
%%% Still not sure exactly what this gives - symmetry and factor for sensitivity by size
simFineTest_freqs = 7.5:0.5:9.5;
%%% Not neccersary to test all powers for 8.5, can combine with previous
simFineTest_powers = round(10.^([1.6 1.7 1.8 2])) % Given threshold at 45
simFineTest_time = simFreqTest_time;

fine_freqCols = winter(length(simFineTest_freqs));
fine_powerCols = cool(length(simFineTest_powers));

% Phase 1.4 - identify time response at low power
%%% Could also do this across combination of frequency and power - not sure of value?
%%% If very short on 1st, probably not neccersary.
simTimeTest_freqs = 8.5;
simTimeTest_powers = round(10.^([1.75 2.25 2.75])); 
%%% In practice, will use previous for 1000
simTimeTest_times = [0.001 0.01 0.1 1 10 100 1000]; 

time_powerCols = cool(length(simTimeTest_powers));
time_timeCols = copper(length(simTimeTest_times));

% To record sensitivity results - should set references dynamically...
sensitivity_FreqRef = 1:0.25:15;
sensitivity_InactRef = 0:5:100;
sensitivity_PowerRef = 1:0.05:4;

inactNoiseThresh = 50;

% Changed all to be indexed by power, inact index doesn't work for indivdiual 
% bulk is basically replicatting line plots
sensitivityBulk = zeros(length(sensitivity_PowerRef), length(sensitivity_FreqRef));
sensitivityIndividual = zeros(length(sensitivity_PowerRef), length(sensitivity_FreqRef), length(influenzaSize));
sensitivityBulk_plaque = zeros(length(sensitivity_PowerRef), length(sensitivity_FreqRef));

%% Make simulated plots
% Set up image
freqRange = min(simFreqTest_freqs):0.2:max(simFreqTest_freqs);
freqRangeFine = 0:0.05:max(simFreqTest_freqs);
powerRange = round(10.^(0:0.05:3)); 0:10:1000;
powerRangeFine = round(10.^(0:0.005:3)); % 0:1:1000; 

markerSize = 6;

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

        if useWeibullPower
            if log10(simFreqTest_power) >= powerWeibullThreshold
                powerVal = 1-exp(-powerWeibullAlpha * ...
                    (log10(simFreqTest_power)-powerWeibullThreshold)^powerWeibullBeta);
            else
               powerVal = 0; 
            end
        else
            powerVal = powerLinearA*(log10(simFreqTest_power))+powerLinearB;
        end
    else
        powerVal = 0;
    end

    % Scale for time, again, should always be 1
    powerScale = log10(simFreqTest_power)/timeWeibullPowerCenter;  
    powerScale = 1 + (powerScale-1)*timeWeibullPowerScale;
    
    if powerScale < 0.01
        powerScale = 0.01;
    end
    
    if log10(simFreqTest_time*powerScale) >= timeWeibullThreshold                   
       timeVal = 1-exp(-timeWeibullAlpha* ...
                (log10(simFreqTest_time*powerScale)-timeWeibullThreshold)^timeWeibullBeta); 
    else
        timeVal = 0;
    end
    
    totalVal = timeVal * powerVal * curveVal;
    totalVal(totalVal > 100) = 100;
    totalVal(totalVal < 0) = 0;
        
    simFreqTest_inactRef(i,:) = totalVal;
    simFreqTest_inact(i,:) = totalVal + absStd*randn(nReps,1);
    
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

%simPowerTest_freqs = freqCurve.b1;



inactRatioBySize = zeros(length(simFreqTest_freqs), nrepsCount, length(influenzaSize));

inactFractionBySize = zeros(length(simFreqTest_freqs), nrepsCount, length(influenzaSize));

activeBySize = zeros(length(simFreqTest_freqs), nrepsCount, length(influenzaSize));

totalInact = zeros(length(simFreqTest_freqs), nrepsCount);

maxInactInd = zeros(length(simFreqTest_freqs), nrepsCount); % if more than 1 equal, will just be one 

maxInactBounds = zeros(length(simFreqTest_freqs), nrepsCount, 2);

% Look at distributions - take for first result
figure; 

refPlotted = zeros(5,1);

%%% Make a function for this, will use later
for i = 1:length(simFreqTest_freqs)
    
    if round(simFreqTest_freqs(i)/4) > 0
        freqPlot = round(simFreqTest_freqs(i)/4);
    else
        freqPlot = 1;
    end
    
    if refPlotted(freqPlot) == 0
        subplot(3,5,5+freqPlot); hold on
        plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-k', 'linewidth', 2)
        refPlotted(freqPlot) = 1;
    end
    
    for j = countVec
        tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

        tempSize_samples = sort(tempSize_samples, 'descend');

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

        samplesIntact = ones(testCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(testCountNum*simFreqTest_inactRef(i,1)/100);
        
        if numToInact > 1
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

        % Strictly speaking, to match equivelent of usual Ct/C0 doesn't require distribution of inactivated
            %%% However, should know proper ratio, which may require fitting sides of active distribution
            zeroRatioInds = find(influenzaSize_dist/initialCountNum*100 < stepPercent);

            zeroFractionInds = find((activeDistDist+inactiveDist)/testCountNum*100 < stepPercent);

            inactRatioBySize(i,j,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

            inactRatioBySize(i,j,zeroRatioInds) = 0;

            inactFractionBySize(i,j,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

            inactFractionBySize(i,j,zeroFractionInds) = 0;

            activeBySize(i,j,:) = activeDistDist/testCountNum*100;

            totalInact(i,j) = (1-sum(activeDistDist)/testCountNum)*100;

            % get inds and ranges on inactivation spectrum
            % find max
            referenceRange = find(inactRatioBySize(i,j,:) == max(inactRatioBySize(i,j,:)));
            
            [maxInact, tempInd] = max(inactRatioBySize(i,j,referenceRange));
            
            maxInactInd(i,j) = referenceRange(tempInd);
            
            tempInactInds = referenceRange(inactRatioBySize(i,j,referenceRange) == maxInact);
            
            maxInactBounds(i,j,:) = influenzaSize(tempInactInds([1 end]));    
            
            
            
            % plotting spectra - counted of inactivated
            subplot(3,5,5+freqPlot); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', freq_freqCols(i,:))

            % inactivation spectrum
            subplot(3,5,2*5+freqPlot); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,:), [4 3 2 1]), 'color', freq_freqCols(i,:));
            ylim([-50 100])
    end
end

subplot(3,5,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
for i = 1:length(simFreqTest_freqs)
    toPlot = setxor(countVec, 1:nReps);
    plot(simFreqTest_FreqRef(i,toPlot), simFreqTest_inact(i,toPlot), ...
        'o', 'markersize', markerSize, 'color', freq_freqCols(i,:));

    plot(simFreqTest_FreqRef(i,countVec), simFreqTest_inact(i,countVec), ...
        '*', 'markersize', markerSize, 'color', freq_freqCols(i,:));
end
xlim([1 20]); ylim([-10 110])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);

plot(simFreqTest_FreqRef(:,1), simFreqTest_inactRef(:,1), 'm-')

% inact from counted intact summed
subplot(3,5,4); hold on
for i = 1:length(simFreqTest_freqs)
    plot(totalInact(i,countVec),  simFreqTest_inact(i,countVec), 'o', 'markersize', markerSize, 'color', freq_freqCols(i,:))
end

line([0 100], [0 100])
ylim([0 100]); xlim([0 100]); 

% Inact from counted both across range
subplot(3,5,2); hold on

for i = 1:length(simFreqTest_freqs)
    plot(simFreqTest_FreqRef(i,countVec), inactFractionBySize(i,countVec, maxInactInd(i,countVec)),...
        'x', 'markersize', markerSize, 'color', freq_freqCols(i,:))
end
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
xlim([1 20])

% Inact from counted intact across range
subplot(3,5,3); hold on
for i = 1:length(simFreqTest_freqs)
        plot(simFreqTest_FreqRef(i,countVec), inactRatioBySize(i,countVec,maxInactInd(i,countVec)),...
            'x', 'markersize', markerSize, 'color', freq_freqCols(i,:))
end
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
xlim([1 20])

plot(influenzSize_resonances, influenzaSize_dist/max(influenzaSize_dist)*100, '-b', 'linewidth', 2)

% Plot inactivation width

%%% Could use thickness or color to indicate number inactivated
subplot(3,5,5); hold on
for i = 1:length(simFreqTest_freqs)
    for j = countVec
        if maxInactBounds(i,j,1) ~= maxInactBounds(i,j,2)
            line(simFreqTest_FreqRef(i,j)*[1 1], permute(maxInactBounds(i,j,:), [3 2 1])*10^9,...
                'color', freq_freqCols(i,:), 'linewidth', 2)
        else
            if inactRatioBySize(i,j,maxInactInd(i,j)) == 100
                plot(simFreqTest_FreqRef(i,j), maxInactBounds(i,j,1)*10^9,...
                    '.', 'markersize', 8, 'color', freq_freqCols(i,:))
            elseif inactRatioBySize(i,j,maxInactInd(i,j)) > inactNoiseThresh
                plot(simFreqTest_FreqRef(i,j), maxInactBounds(i,j,1)*10^9,...
                    'o', 'markersize', 4, 'color', freq_freqCols(i,:))
            end
        end
    end
end
ylim([80 120])
xlim([1 20])

freqInactRatioBySize = inactRatioBySize;



%%% Store sensitivity and plot
    % Still neccersary?
for i = 1:length(simFreqTest_freqs)
    if simFreqTest_freqs(i) > min(sensitivity_FreqRef) & simFreqTest_freqs(i) < max(sensitivity_FreqRef)
        [~, freqInd] = min(abs(sensitivity_FreqRef - simFreqTest_freqs(i)));

        [~, powerInd] = min(abs(sensitivity_PowerRef - log10(simFreqTest_power)));

        sensitivityBulk(powerInd, freqInd) = mean(totalInact(i,countVec));

        if mean(simFreqTest_inact(i,:)) > 0
            sensitivityBulk_plaque(powerInd, freqInd) = mean(simFreqTest_inact(i,:));
        end
        
        for j = 1:length(influenzaSize)
            if mean(inactRatioBySize(i,countVec,j)) > inactNoiseThresh
                sensitivityIndividual(powerInd, freqInd, j) = mean(inactRatioBySize(i,countVec,j));
            end
        end
    end
end

sensitivityFig = figure;

subplot(4,5,1); hold on
imshow(log10(sensitivityBulk)/2)

subplot(4,5,2); hold on
imshow(log10(sensitivityBulk_plaque)/2)

sizeToPlot = [90 100 110];
for i = 1:length(sizeToPlot)
    subplot(4,5,2+i); hold on
    
    sizeInd = find(influenzaSize*10^9 == sizeToPlot(i));
    
    imshow(log10(sensitivityIndividual(:,:,sizeInd))/2)
end

%%% Removed image code

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

        % Scale for power
        if simPowerTest_powers(j) > powerThresholdInterp_powerTest(i)
            if useWeibullPower
                if log10(simPowerTest_powers(j)) >= powerWeibullThreshold
                    powerVal = 1-exp(-powerWeibullAlpha * ...
                        (log10(simPowerTest_powers(j))-powerWeibullThreshold)^powerWeibullBeta);
                else
                   powerVal = 0; 
                end
            else
                powerVal = powerLinearA*(log10(simPowerTest_powers(j)))+powerLinearB;    
            end
        else
            powerVal = 0;
        end
        
        % Scale for time
        powerScale = log10(simPowerTest_powers(j))/timeWeibullPowerCenter;   
        powerScale = 1 + (powerScale-1)*timeWeibullPowerScale;

        if powerScale < 0.01
            powerScale = 0.01;
        end
        
        if log10(simPowerTest_time*powerScale) >= timeWeibullThreshold                   
           timeVal = 1-exp(-timeWeibullAlpha* ...
                    (log10(simPowerTest_time*powerScale)-timeWeibullThreshold)^timeWeibullBeta); 
        else
            timeVal = 0;
        end

        totalVal = timeVal * powerVal * curveVal;
        totalVal(totalVal > 100) = 100;     
        totalVal(totalVal < 0) = 0;

        simPowerTest_inactRef(i,j,:) = totalVal;
        simPowerTest_inact(i,j,:) = totalVal + absStd*randn(nReps,1);
        
        simPowerTest_freqRef(i,j,:) = simPowerTest_freqs(i);
        simPowerTest_powerRef(i,j,:) = simPowerTest_powers(j);
    end
end

simPowerTest_inact(simPowerTest_inact > 100) = 100;

% Plot power curve
figure; 
subplot(4,length(simPowerTest_powers),1); hold on

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        toPlot = setxor(countVec, 1:nReps);

        plot3(log10(permute(simPowerTest_powerRef(i,j,toPlot), [3 2 1])), permute(simPowerTest_inact(i,j,toPlot), [3 2 1]), permute(simPowerTest_freqRef(i,j,toPlot), [3 2 1]), ...
            'o', 'markersize', markerSize, 'color', power_freqCols(i,:))

        plot3(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(simPowerTest_inact(i,j,countVec), [3 2 1]), permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), ...
            '*', 'markersize', markerSize, 'color', power_freqCols(i,:))
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

plot(log10(permute(simPowerTest_powerRef(1,:,1), [3 2 1])), permute(simPowerTest_inactRef(1,:,1), [3 2 1]), 'm', 'linewidth', 2)

% Plot freq curve
subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+1); hold on

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        toPlot = setxor(countVec, 1:nReps); 

        plot3(permute(simPowerTest_freqRef(i,j,toPlot), [3 2 1]), permute(simPowerTest_inact(i,j,toPlot), [3 2 1]), log10(permute(simPowerTest_powerRef(i,j,toPlot), [3 2 1])), ...
            'o', 'markersize', markerSize, 'color', power_powerCols(j,:))

        plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(simPowerTest_inact(i,j,countVec), [3 2 1]), log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), ...
            '*', 'markersize', markerSize, 'color', power_powerCols(j,:))
    end
end

plot(simPowerTest_freqRef(:,1,1), simPowerTest_inactRef(:,1,1), 'm', 'linewidth', 2)

plot(simPowerTest_freqRef(:,end,1), simPowerTest_inactRef(:,end,1), 'm', 'linewidth', 2)



inactRatioBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount, length(influenzaSize));

inactFractionBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount, length(influenzaSize));

activeBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount, length(influenzaSize));

totalInact = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount);


inactInds = cell(length(simPowerTest_freqs), length(simPowerTest_powers));

maxInactInd = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount); % if more than 1 equal, will just be one 

maxInactBounds = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount, 2);


%%% Main to do here
%%%     1. scale active virus particle dist sides to match original dist
%%%     2. Pick drop in dist/peak innact and track as it gets lower - otherwise centre point needs to be regularly adjusted
            % May move by dist

%%% Make a function for this, will use later
for i = 1:length(simPowerTest_freqs)
    for j = fliplr(1:length(simPowerTest_powers))

        subplot(4,length(simPowerTest_powers),2*length(simPowerTest_powers)+j); hold on
        plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-k', 'linewidth', 2)

        for k = countVec
            tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

            tempSize_samples = sort(tempSize_samples, 'descend');

            tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

            samplesIntact = ones(testCountNum, 1);

            % Start with reference before we figure out noise... 
            numToInact = round(testCountNum*simPowerTest_inactRef(i,j,k)/100);

            if numToInact > 1 
                % Find nearest in sample
                [~, minInd] = min(abs(tempSize_freqs-simPowerTest_freqs(i)));

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

                        [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-simPowerTest_freqs(i)));

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

            inactRatioBySize(i,j,k,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

            inactRatioBySize(i,j,k,zeroRatioInds) = 0;

            inactFractionBySize(i,j,k,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

            inactFractionBySize(i,j,k,zeroFractionInds) = 0;

            activeBySize(i,j,k,:) = activeDistDist/testCountNum*100;

            totalInact(i,j,k) = (1-sum(activeDistDist)/testCountNum)*100;

            % get inds and ranges on inactivation spectrum
            % find max
            if length(simPowerTest_powers) == j
                referenceRange = find(inactRatioBySize(i,j,k,:) == max(inactRatioBySize(i,j,k,:)));
                
                referenceRange = union(inactInds{i,j}, referenceRange);
            else
                referenceRange = inactInds{i, j+1};
            end
            
            [maxInact, tempInd] = max(inactRatioBySize(i,j,k,referenceRange));
            
            maxInactInd(i,j,k) = referenceRange(tempInd);
            
            % Note that union is used so that this would grow across multiple count replicats
            tempInactInds = referenceRange(inactRatioBySize(i,j,k,referenceRange) == maxInact);
            
            inactInds{i,j} = union(tempInactInds, inactInds{i,j});
            
            maxInactBounds(i,j,k,:) = influenzaSize(tempInactInds([1 end]));    
            
            
            
            % plotting spectra - counted of inactivated
            subplot(4,length(simPowerTest_powers),2*length(simPowerTest_powers)+j); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', power_freqCols(i,:))

            % inactivation spectrum
            subplot(4,length(simPowerTest_powers),3*length(simPowerTest_powers)+j); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,:), [4 3 2 1]), 'color', power_freqCols(i,:))
            ylim([-50 100])
        end
    end
end

%%% Plot factored by power
    % inact from counted intact summed
    subplot(4,length(simPowerTest_powers),4); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            plot(permute(totalInact(i,j,countVec), [4 3 2 1]),  permute(simPowerTest_inact(i,j,countVec), [3 2 1]), ...
                'o', 'markersize', markerSize, 'color', power_freqCols(i,:))
        end
    end

    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]); 


    % Inact from counted both across range
    subplot(4,length(simPowerTest_powers),2); hold on

    for i = 1:length(simPowerTest_freqs) 
        for j = 1:length(simPowerTest_powers)
            plot3(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4, 3 2 1]), permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), ...
                'x', 'markersize', markerSize, 'color', power_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Inact from counted intact across range
    subplot(4,length(simPowerTest_powers),3); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            plot3(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]), permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), ...
                'x', 'markersize', markerSize, 'color', power_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simPowerTest_powers),5); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            for k = countVec
                if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2)
                    line(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9, permute(simPowerTest_freqRef(i,j,countVec), [3 2 1])*[1 1], ...
                        'color', power_freqCols(i,:), 'linewidth', 2)
                    else
                        if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100
                            plot3(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), ...
                                '.', 'markersize', 8, 'color', power_freqCols(i,:))
                        elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) > inactNoiseThresh
                            plot3(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), ...
                                'o', 'markersize', 4, 'color', power_freqCols(i,:))
                        end
                    end
            end
        end
    end
    ylim([80 120])
    xlim([1 3])
    
%%% Plot factored by freq
    % inact from counted intact summed
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+4); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            plot(permute(totalInact(i,j,countVec), [4 3 2 1]),  permute(simPowerTest_inact(i,j,countVec), [3 2 1]), ...
                'o', 'markersize', markerSize, 'color', power_powerCols(j,:))
        end
    end

    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]); 


    % Inact from counted both across range
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+2); hold on

    for i = 1:length(simPowerTest_freqs) 
        for j = 1:length(simPowerTest_powers) 
            plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4, 3 2 1]), log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), ...
                'x', 'markersize', markerSize, 'color', power_powerCols(j,:))
        end
    end
    
    % Inact from counted intact across range
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+3); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers) 
            plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]), log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), ...
                'x', 'markersize', markerSize, 'color', power_powerCols(j,:))
        end
    end

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+5); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            for k = countVec
                if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2) 
                    line(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1])*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9, log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1]))*[1 1], ...
                        'color', power_powerCols(j,:), 'linewidth', 2)
                else
                    if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100  
                        plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), ...
                            '.', 'markersize', 8, 'color', power_powerCols(j,:))
                    elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) > inactNoiseThresh
                        plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), ...
                            'o', 'markersize', 4, 'color', power_powerCols(j,:))
                    end
                end
            end
        end
    end
    ylim([80 120])
    
powerInactRatioBySize = inactRatioBySize;
    
%%% Store sensitivity and plot
for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        if simPowerTest_freqs(i) > min(sensitivity_FreqRef) & simPowerTest_freqs(i) < max(sensitivity_FreqRef)
            [~, freqInd] = min(abs(sensitivity_FreqRef - simPowerTest_freqs(i)));

            [~, powerInd] = min(abs(sensitivity_PowerRef - log10(simPowerTest_powers(j))));

            if sensitivityBulk(powerInd, freqInd) == 0
                sensitivityBulk(powerInd, freqInd) = mean(totalInact(i,j,countVec));
            else
                sensitivityBulk(powerInd, freqInd) = min([mean(totalInact(i,j,countVec)) sensitivityBulk(powerInd, freqInd)]);
            end
            
            if mean(simPowerTest_inact(i,j,:)) > 0
                if sensitivityBulk_plaque(powerInd, freqInd) == 0
                    sensitivityBulk_plaque(powerInd, freqInd) = mean(simPowerTest_inact(i,j,:));
                else
                    sensitivityBulk_plaque(powerInd, freqInd) = min([mean(simPowerTest_inact(i,j,:)) sensitivityBulk_plaque(powerInd, freqInd)]);
                end
            end

            for k = 1:length(influenzaSize)
                if mean(inactRatioBySize(i,j,countVec,k)) > inactNoiseThresh
                    if sensitivityIndividual(powerInd, freqInd, k) == 0
                        sensitivityIndividual(powerInd, freqInd, k) = mean(inactRatioBySize(i,j,countVec,k));
                    else
                        sensitivityIndividual(powerInd, freqInd, k) = min([mean(inactRatioBySize(i,j,countVec,k)) sensitivityIndividual(powerInd, freqInd, k)]);
                    end
                end
            end
        end
    end
end

figure(sensitivityFig);

subplot(4,5,5+1); hold on
imshow(log10(sensitivityBulk)/2)

subplot(4,5,5+2); hold on
imshow(log10(sensitivityBulk_plaque)/2)

sizeToPlot = [90 100 110];
for i = 1:length(sizeToPlot)
    subplot(4,5,5+2+i); hold on
    
    sizeInd = find(influenzaSize*10^9 == sizeToPlot(i));
    
    imshow(log10(sensitivityIndividual(:,:,sizeInd))/2)
end


inds = find(sensitivityIndividual);
[x, y, z] = ind2sub(size(sensitivityIndividual), inds);

figure; subplot(1,2,1)
plot3(sensitivity_PowerRef(x), sensitivity_FreqRef(y), influenzaSize(z)*10^9, '.')

subplot(1,2,2)
plot3(sensitivity_PowerRef(x), sensitivity_FreqRef(y), influenzSize_resonances(z), '.')

%% Solver using both freq and power scan

removeUnused = 1;

swapInReference = 1;

predictionFreqs = sort([simFreqTest_freqs simPowerTest_freqs]);
predictionPowers = simPowerTest_powers;

predictionSizes = influenzaSize;
predictionSizes_dist = influenzaSize_dist/sum(influenzaSize_dist)*100;

predictionTime = simPowerTest_time;

% Get matched power across freq
if length(powerThreshold) == 1
    powerThresholdInterp_pred = powerThreshold*ones(length(predictionFreqs),1);
else
    % Interp if more than one
    powerThresholdInterp_pred = interp1(powerThresholdFreqs, powerThreshold, predictionFreqs, 'linear',0);
    powerThresholdInterp_pred(powerThresholdInterp_pred == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        predictionFreqs(powerThresholdInterp_pred == 0), 'nearest','extrap');
end

predictedInactivation = zeros(length(predictionFreqs),length(predictionPowers));

%%% Testing using known curves, need to change to fitted.
for i = 1:length(predictionFreqs)
    curveVal = curveMax*exp(-(predictionFreqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(predictionPowers)
        if predictionPowers(j) > powerThresholdInterp_pred(i)
            if useWeibullPower
                if log10(predictionPowers(j)) >= powerWeibullThreshold
                    powerVal = 1-exp(-powerWeibullAlpha * ...
                        (log10(predictionPowers(j))-powerWeibullThreshold)^powerWeibullBeta);
                else
                   powerVal = 0; 
                end
            else
                powerVal = powerLinearA*(log10(predictionPowers(j)))+powerLinearB;    
            end
        else
            powerVal = 0;
        end

        % Scale for time
        powerScale = log10(predictionPowers(j))/timeWeibullPowerCenter;   
        powerScale = 1 + (powerScale-1)*timeWeibullPowerScale;

        if powerScale < 0.01
            powerScale = 0.01;
        end

        if log10(predictionTime*powerScale) >= timeWeibullThreshold                   
           timeVal = 1-exp(-timeWeibullAlpha* ...
                    (log10(predictionTime*powerScale)-timeWeibullThreshold)^timeWeibullBeta); 
        else
            timeVal = 0;
        end

        totalVal = timeVal * powerVal * curveVal;
        totalVal(totalVal > 100) = 100;     
        totalVal(totalVal < 0) = 0;

        predictedInactivation(i,j) = totalVal;
    end
end

% Remove low points from power
if removeUnused
    predictInactThresh = 0.1; %2;

    toRemove = zeros(length(predictionPowers),1,'logical');

    for i = 1:length(predictionPowers)
        if all( predictedInactivation(:,i) < predictInactThresh)
            toRemove(i) = 1;
        end
    end

    powersRemoved = sum(toRemove);
    
    predictionPowers(toRemove) = [];
    predictedInactivation(:,toRemove) = [];
end

thresholdArray = zeros(length(predictionFreqs), length(influenzaSize));

% initalize threshold array with already measured values

for i = 1:length(simFreqTest_freqs)
    freqInd = find(simFreqTest_freqs(i) == predictionFreqs);
    
    if ~isempty(freqInd)
        
        meanInacts = mean(freqInactRatioBySize(i,countVec,:),2);
        
        sizeInds = find(meanInacts > inactNoiseThresh);
        
        [inds] = sub2ind(size(thresholdArray), freqInd*ones(length(sizeInds),1), sizeInds);
        
%         thresholdArray(inds) = simFreqTest_power;
        thresholdArray(inds) = length(simPowerTest_powers)-powersRemoved;
    end
end

for i = 1:length(simPowerTest_freqs)
    freqInd = find(simPowerTest_freqs(i) == predictionFreqs);
    
    if ~isempty(freqInd)
        for j = fliplr(1:length(simPowerTest_powers))
            
            meanInacts = mean(powerInactRatioBySize(i,j,countVec,:),2);

            sizeInds = find(meanInacts > inactNoiseThresh);

            [inds] = sub2ind(size(thresholdArray), freqInd*ones(length(sizeInds),1), sizeInds);

%             thresholdArray(inds) = simPowerTest_powers(j);
            thresholdArray(inds) = j-powersRemoved;
        end
    end
end

% Get reference inactivation for unused
sizesToUse = 1:length(influenzaSize);    

inactRatioBySize_reference = zeros(length(predictionFreqs), length(predictionPowers), length(predictionSizes));

totalInact = zeros(length(predictionFreqs), length(predictionPowers));

lostInact = zeros(length(predictionFreqs), length(predictionPowers));
keptInact = zeros(length(predictionFreqs), length(predictionPowers));

for i = 1:length(predictionFreqs)
    for j = fliplr(1:length(predictionPowers))
        tempSize_samples = influenzaSize_samples;

        tempSize_samples = sort(tempSize_samples, 'descend');

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

        samplesIntact = ones(testCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(testCountNum*predictedInactivation(i,j)/100);

        if numToInact > 1 
            % Find nearest in sample
            [~, minInd] = min(abs(tempSize_freqs-predictionFreqs(i)));

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

                    [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-predictionFreqs(i)));

                    minInd = indRef(tempInd);
                end

                samplesIntact(minInd) = 0;

                numInactivated = numInactivated + 1;
            end
        end

        % Get hist of both
        activeDistDist = hist(tempSize_samples(samplesIntact == 1), influenzaSize);

        % Strictly speaking, to match equivelent of usual Ct/C0 doesn't require distribution of inactivated
        %%% However, should know proper ratio, which may require fitting sides of active distribution

        inactRatioBySize_reference(i,j,:) = (1 - (activeDistDist(sizesToUse)/predictionCountNum)./...
            (influenzaSize_dist(sizesToUse)/predictionCountNum))*100;
        
        % for testing
        totalInact(i,j) = (1 - sum(activeDistDist)/initialCountNum)*100;
        
        %%% Need to store data to calculate later... 
        lostInact(i,j) = (1 - sum(activeDistDist(sizesRemoved))/ignoredCountNum)*100;
        
        keptInact(i,j) = (1 - sum(activeDistDist(sizesToUse))/predictionCountNum)*100;  
    end
end

lostInact(isnan(lostInact)) = 0;

% Collapse to threshold map

thresholdArray_reference = zeros(length(predictionFreqs), length(predictionSizes));

for i = 1:length(predictionFreqs)
    for j = fliplr(1:length(predictionPowers))

        sizeInds = find(inactRatioBySize_reference(i,j,:) > 50); %inactNoiseThresh

        [inds] = sub2ind(size(thresholdArray_reference), i*ones(length(sizeInds),1), sizeInds);

        thresholdArray_reference(inds) = j;
    end
end

if swapInReference
    thresholdArray = thresholdArray_reference; 
    
    predictedInactivation = totalInact; 
end

% Remove from sizes and freqs
toRemove = zeros(length(predictionSizes),1,'logical');

predictionSizes_distTemp = predictionSizes_dist;

if removeUnused
    toRemoveTemp = zeros(length(predictionFreqs),1, 'logical');

    for i = 1:length(predictionFreqs)
        if all( thresholdArray(i,:) == 0)
            toRemoveTemp(i) = 1;
        end
    end

    predictionFreqs(toRemoveTemp) = [];
    predictedInactivation(toRemoveTemp,:) = [];
    thresholdArray(toRemoveTemp,:) = [];
    
    thresholdArray_reference(toRemoveTemp,:) = [];
    inactRatioBySize_reference(toRemoveTemp,:) = [];
    totalInact(toRemoveTemp,:) = [];
    lostInact(toRemoveTemp,:) = [];
    keptInact(toRemoveTemp,:) = [];
    
    for i = 1:length(predictionSizes)
        if all( thresholdArray(:,i) == 0)
            toRemove(i) = 1;
        end
    end

    sizesToUse(toRemove) = [];
    predictionSizes(toRemove) = [];
    predictionSizes_dist(toRemove) = [];
    thresholdArray(:,toRemove) = [];
    
    thresholdArray_reference(:,toRemove) = [];
    inactRatioBySize_reference(:,toRemove) = [];
end

sizesRemoved = find(toRemove);
predictionCountNum = sum(influenzaSize_dist(sizesToUse));
ignoredCountNum = sum(influenzaSize_dist(toRemove));

if any(influenzaSize_dist(sizesToUse)/predictionCountNum*100 < stepPercent)
    warning('small sizes are present')
end

% Lock in lines that are already measured
tempArray = thresholdArray;
for i = 1:length(simPowerTest_freqs)
    freqInd = find(simPowerTest_freqs(i) == predictionFreqs);
    tempArray(freqInd,:) = 0;
end    
    
pointsToSolve = find(tempArray == length(simPowerTest_powers)-powersRemoved);

% testing - dif predicted to total - should be similar
figure; 
subplot(5,3,1); imshow(predictedInactivation/100)

subplot(5,3,2); imshow(totalInact/100)

subplot(5,3,3); imshow(abs(predictedInactivation-totalInact)*10)
        
% testing - dif predicted to sum of kept and lost inact - should be similar
% similar
testInactivation = lostInact*ignoredCountNum/initialCountNum + keptInact*predictionCountNum/initialCountNum;

subplot(5,3,4); imshow(abs(totalInact-testInactivation)*100)

subplot(5,3,5); imshow(testInactivation/100)

subplot(5,3,6); imshow(abs(predictedInactivation-testInactivation)*10)

% testing - dif predictied to kept inact - will be dif

subplot(5,3,8); imshow(keptInact/100)

subplot(5,3,9); imshow(abs(predictedInactivation-keptInact)*10)

% testing - dif predictied to sum of kept inact across size - will be dif

testInactivation = size(predictedInactivation);

for i = 1:length(predictionFreqs)
    for j = 1:length(predictionPowers)

        testInactivation(i,j) = sum(permute(inactRatioBySize_reference(i,j,:), [2 3 1])/100 .* predictionSizes_dist); %...
%             influenzaSize_dist(sizesToUse)/predictionCountNum*100);
        
    end
end

subplot(5,3,10); imshow(abs(testInactivation - keptInact)*100)

subplot(5,3,11); imshow(testInactivation/100)

subplot(5,3,12); imshow(abs(predictedInactivation-testInactivation)*10)

% testing - scaled dif to kept inact

subplot(5,3,14); imshow(abs(predictedInactivation.*predictionCountNum/initialCountNum-keptInact)*10)

subplot(5,3,15); imshow(abs(predictedInactivation-keptInact./(predictionCountNum/initialCountNum))*10)

tempPoints = find(thresholdArray_reference);

% testing - fn
fun = @(x)inactivationError(x, (thresholdArray), predictionSizes_dist, ... 
    keptInact, tempPoints, 1:length(predictionPowers));

[sse, ~, fnInact] = fun(thresholdArray_reference(tempPoints));

testInactivation2 = size(predictedInactivation);

for i = 1:length(predictionFreqs)
    for j = 1:length(predictionPowers)

        inds = find(permute(inactRatioBySize_reference(i,j,:), [3 2 1]) > 50);
        
        testInactivation2(i,j) = sum(predictionSizes_dist(inds));
    end
end

%should match
[sse sum((predictedInactivation(:)-totalInact(:)).^2)]

figure;
subplot(1,2,1)
imshow(abs(fnInact - testInactivation)/10)

subplot(1,2,2)
imshow(abs(fnInact - testInactivation2)*10)

%%% Solver

% fun = @(x)inactivationError(x, log10(thresholdArray), predictionSizes_dist, ... 
%     predictedInactivation, pointsToSolve, log10(predictionPowers)); 

fun = @(x)inactivationError(x, (thresholdArray), predictionSizes_dist, ... 
    predictedInactivation, pointsToSolve, 1:length(predictionPowers));

allowBroad = 1;

con = @(x)inactivationConstraints(x, (thresholdArray), pointsToSolve, allowBroad); %log10

% startVals = log10(max(simPowerTest_powers)*ones(length(pointsToSolve),1)); %log10
% ubd = log10(max(simPowerTest_powers)*ones(length(pointsToSolve),1)); %log10
% lbd = log10(min(simPowerTest_powers)*ones(length(pointsToSolve),1)); %log10

startVals = (length(predictionPowers)*ones(length(pointsToSolve),1)); 
ubd = (length(predictionPowers)*ones(length(pointsToSolve),1)); 
lbd = (ones(length(pointsToSolve),1)); 

% Seems to work quite well with integer points on pattern search
    % Recommend is to us ga, but doesn't allow equality constraints
    
%%% Set up flag to switch between integer or log values     
    
% can set to use parallel
patSearOpts = optimoptions('patternsearch','Display','iter', 'FunctionTolerance', 1e-6, ...
            'MaxFunctionEvaluations', 10^6, 'MaxIterations', 10^3, 'ScaleMesh', 'off', 'TolMesh', 0.9);
    
vals = patternsearch(fun, startVals, [], [], [], [], lbd, ubd, con, patSearOpts);

% Never seems to work with fmincon

[a, b, c] = fun(vals);

a

[c, ceq] = con(vals);

% [(1:length(ceq))' 10.^b' ceq]'

b(b > 0) = predictionPowers(ceil(b(b > 0)));

[(1:size(b,2))' b' ceq(1:size(b,2)) ceq(size(b,2)+1:end)]'
% [(1:length(ceq))' thresholdArray' ceq]'

thresholdTemp = thresholdArray;
thresholdTemp(pointsToSolve) = vals;

figure; 
subplot(1,3,1); imshow(thresholdArray_reference/length(predictionPowers))

subplot(1,3,2); imshow(thresholdTemp/length(predictionPowers))

subplot(1,3,3); imshow(abs(thresholdArray_reference-thresholdTemp)/length(predictionPowers))

%%% Want to get confidence intervals for each fit

%% Phase 1.3 scan across freq

simFineTest_inact = zeros(length(simFineTest_freqs), length(simFineTest_powers), nReps);
simFineTest_inactRef = zeros(length(simFineTest_freqs), length(simFineTest_powers), nReps);
simFineTest_freqRef = zeros(length(simFineTest_freqs), length(simFineTest_powers), nReps);
simFineTest_powerRef = zeros(length(simFineTest_freqs), length(simFineTest_powers), nReps);

% Get matched roots across freq
if length(powerThreshold) == 1
    powerThresholdInterp_powerTest = powerThreshold*ones(length(simFineTest_freqs),1);
else
    powerThresholdInterp_powerTest = interp1(powerThresholdFreqs, powerThreshold, simFineTest_freqs, 'linear',0);
    powerThresholdInterp_powerTest(powerThresholdInterp_powerTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        simFineTest_freqs(powerThresholdInterp_powerTest == 0), 'nearest','extrap');
end

for i = 1:length(simFineTest_freqs)

    curveVal = curveMax*exp(-(simFineTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simFineTest_powers)

        % scale for power
        if simFineTest_powers(j) > powerThresholdInterp_powerTest(i)
            if useWeibullPower
                if log10(simFineTest_powers(j)) >= powerWeibullThreshold
                    powerVal = 1-exp(-powerWeibullAlpha * ...
                        (log10(simFineTest_powers(j))-powerWeibullThreshold)^powerWeibullBeta);
                else
                   powerVal = 0; 
                end
            else
                powerVal = powerLinearA*(log10(simFineTest_powers(j)))+powerLinearB;    
            end
        else
            powerVal = 0;
        end
        
        % Scale for time
        powerScale = log10(simFineTest_powers(j))/timeWeibullPowerCenter;   
        powerScale = 1 + (powerScale-1)*timeWeibullPowerScale;

        if powerScale < 0.01
            powerScale = 0.01;
        end
        
        if log10(simFineTest_time*powerScale) >= timeWeibullThreshold                   
           timeVal = 1-exp(-timeWeibullAlpha* ...
                    (log10(simFineTest_time*powerScale)-timeWeibullThreshold)^timeWeibullBeta); 
        else
            timeVal = 0;
        end

        totalVal = timeVal * powerVal * curveVal;
        totalVal(totalVal > 100) = 100;     
        totalVal(totalVal < 0) = 0;

        simFineTest_inactRef(i,j,:) = totalVal;
        simFineTest_inact(i,j,:) = totalVal + absStd*randn(nReps,1);

        simFineTest_freqRef(i,j,:) = simFineTest_freqs(i);
        simFineTest_powerRef(i,j,:) = simFineTest_powers(j);
    end
end

simFineTest_inact(simFineTest_inact > 100) = 100;

% Plot power curve
figure; 
subplot(4,length(simFineTest_freqs),1); hold on

for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)
        toPlot = setxor(countVec, 1:nReps);
        plot3(log10(permute(simFineTest_powerRef(i,j,toPlot), [3 2 1])), permute(simFineTest_inact(i,j,toPlot), [3 2 1]), permute(simFineTest_freqRef(i,j,toPlot), [3 2 1]), ...
            'o', 'markersize', markerSize, 'color', fine_freqCols(i,:));
        
        plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(simFineTest_inact(i,j,countVec), [3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
            '*', 'markersize', markerSize, 'color', fine_freqCols(i,:));
    end
end

plot3(log10(permute(simFineTest_powerRef(end,:,1), [2 1 3])), permute(simFineTest_inactRef(end,:,1), [2 1 3]), permute(simFineTest_freqRef(end,:,1), [2 1 3]),...
    'm', 'linewidth', 2)

plot3(log10(permute(simFineTest_powerRef(1,:,1), [2 1 3])), permute(simFineTest_inactRef(1,:,1), [2 1 3]), permute(simFineTest_freqRef(1,:,1), [2 1 3]),...
    ':m', 'linewidth', 2)

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

%plot time curve
subplot(4,length(simFineTest_freqs),length(simFineTest_freqs)+1); hold on

for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)
            toPlot = setxor(countVec, 1:nReps);
            plot3(permute(simFineTest_freqRef(i,j,toPlot), [3 2 1]), permute(simFineTest_inact(i,j,toPlot), [3 2 1]), log10(permute(simFineTest_powerRef(i,j,toPlot), [3 2 1])), ...
                'o', 'markersize', markerSize, 'color', fine_powerCols(j,:));

            plot3(permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), permute(simFineTest_inact(i,j,countVec), [3 2 1]), log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), ...
                '*', 'markersize', markerSize, 'color', fine_powerCols(j,:));
    end
end

plot3(permute(simFineTest_freqRef(:,end,1), [1 2 3]), permute(simFineTest_inactRef(:,end,1), [1 2 3]), log10(permute(simFineTest_powerRef(:,end,1), [1 2 3])),...
    'm', 'linewidth', 2)

plot3(permute(simFineTest_freqRef(:,1,1), [1 2 3]), permute(simFineTest_inactRef(:,1,1), [1 2 3]), log10(permute(simFineTest_powerRef(:,1,1), [1 2 3])),...
    ':m', 'linewidth', 2)

xlim([ 7 10])

inactRatioBySize = zeros(length(simFineTest_freqs), length(simFineTest_powers), nrepsCount, length(influenzaSize));

inactFractionBySize = zeros(length(simFineTest_freqs), length(simFineTest_powers), nrepsCount, length(influenzaSize));

activeBySize = zeros(length(simFineTest_freqs), length(simFineTest_powers), nrepsCount, length(influenzaSize));

totalInact = zeros(length(simFineTest_freqs), length(simPowerTest_powers), nrepsCount);


inactInds = cell(length(simFineTest_freqs), length(simPowerTest_powers));

maxInactInd = zeros(length(simFineTest_freqs), length(simPowerTest_powers), nrepsCount); % if more than 1 equal, will just be one 

maxInactBounds = zeros(length(simFineTest_freqs), length(simPowerTest_powers), nrepsCount, 2);

%%% Main to do here
%%%     1. scale active virus particle dist sides to match original dist
%%%     2. Pick drop in dist/peak innact and track as it gets lower - otherwise centre point needs to be regularly adjusted
            % May move by dist

%%% Make a function for this, will use later
for i = 1:length(simFineTest_freqs)
    
    subplot(4,length(simFineTest_freqs),2*length(simFineTest_freqs)+i); hold on
    plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-b', 'linewidth', 2)
        
    for j = fliplr(1:length(simFineTest_powers))
        for k = countVec
            tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

            tempSize_samples = sort(tempSize_samples, 'descend');

            tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

            samplesIntact = ones(testCountNum, 1);

            % Start with reference before we figure out noise... 
            numToInact = round(testCountNum*simFineTest_inactRef(i,j,1)/100);

            if numToInact > 1
                % Find nearest in sample
                [~, minInd] = min(abs(tempSize_freqs-simFineTest_freqs(i)));

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

                        [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-simFineTest_freqs(i)));

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

            inactRatioBySize(i,j,k,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

            inactRatioBySize(i,j,k,zeroRatioInds) = 0;

            inactFractionBySize(i,j,k,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

            inactFractionBySize(i,j,k,zeroFractionInds) = 0;

            activeBySize(i,j,k,:) = activeDistDist/testCountNum*100;

            totalInact(i,j,k) = (1-sum(activeDistDist)/testCountNum)*100;
            
             % get inds and ranges on inactivation spectrum
            
            % find max
            if length(simFineTest_powers) == j
                referenceRange = find(inactRatioBySize(i,j,k,:) == max(inactRatioBySize(i,j,k,:)));
                
                referenceRange = union(inactInds{i,j}, referenceRange);
            else
                referenceRange = inactInds{i, j+1};
            end
            
            [maxInact, tempInd] = max(inactRatioBySize(i,j,k,referenceRange));
            
            maxInactInd(i,j,k) = referenceRange(tempInd);
            
            % Note that union is used so that this would grow across multiple count replicats
            tempInactInds = referenceRange(inactRatioBySize(i,j,k,referenceRange) == maxInact);
            
            inactInds{i,j} = union(tempInactInds, inactInds{i,j});
            
            maxInactBounds(i,j,k,:) = influenzaSize(tempInactInds([1 end]));
            
            % plotting spectra - counted of inactivated
            subplot(4,length(simFineTest_freqs),2*length(simFineTest_freqs)+i); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', fine_powerCols(j,:))

            % inactivation spectrum
            subplot(4,length(simFineTest_freqs),3*length(simFineTest_freqs)+i); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,:), [4 3 2 1]), 'color', fine_powerCols(j,:));
            ylim([-50 100])
        end
    end
end

%%% Factored by power
    % inact from counted intact summed
    subplot(4,length(simFineTest_freqs),4); hold on
    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)

            plot(permute(totalInact(i,j,countVec), [4 3 2 1]), permute(simFineTest_inact(i,j,countVec), [3 2 1]), 'o', 'markersize', markerSize, 'color', fine_freqCols(i,:))
        end
    end
    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]);

    % Inact from counted both across range
    subplot(4,length(simFineTest_freqs),2); hold on

    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4 3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
                'x', 'markersize', markerSize, 'color', fine_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Inact from counted intact across range
    subplot(4,length(simFineTest_freqs),3); hold on

    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]), simFineTest_freqRef(i,j,countVec), ...
                'x', 'markersize', markerSize, 'color', fine_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simFineTest_freqs),5); hold on
    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            for k = countVec
                if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2)
                    line(log10(permute(simFineTest_powerRef(i,j,k), [3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9, permute(simFineTest_freqRef(i,j,countVec), [3 2 1])*[1 1],...
                        'color', fine_freqCols(i,:), 'linewidth', 2)
                else
                    if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100
                        plot3(log10(permute(simFineTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, permute(simFineTest_freqRef(i,j,countVec), [3 2 1]),...
                            '.', 'color', fine_freqCols(i,:), 'markersize', 8)
                    elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) > inactNoiseThresh
                        plot3(log10(permute(simFineTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, permute(simFineTest_freqRef(i,j,countVec), [3 2 1]),...
                            'o', 'color', fine_freqCols(i,:), 'markersize', 4)
                    end
                end
            end
        end
    end

    ylim([80 120])
    xlim([1 3])
    
%%% Factored by frequency
    % inact from counted intact summed
    subplot(4,length(simFineTest_freqs),length(simFineTest_freqs)+4); hold on
    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)

            plot(permute(totalInact(i,j,countVec), [4 3 2 1]), permute(simFineTest_inact(i,j,countVec), [3 2 1]), 'o', 'markersize', markerSize, 'color', fine_powerCols(j,:))
        end
    end
    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]);

    % Inact from counted both across range
    subplot(4,length(simFineTest_freqs),length(simFineTest_freqs)+2); hold on

    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            plot3( permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4 3 2 1]), log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), ...
                'x', 'markersize', markerSize, 'color', fine_powerCols(j,:))
        end
    end
    xlim([7 10])

    % Inact from counted intact across range
    subplot(4,length(simFineTest_freqs),length(simFineTest_freqs)+3); hold on

    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            plot3(simFineTest_freqRef(i,j,countVec), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]), log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), ...
                'x', 'markersize', markerSize, 'color', fine_powerCols(j,:)) 
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([7 10])

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simFineTest_freqs),length(simFineTest_freqs)+5); hold on
    for i = 1:length(simFineTest_freqs)
        for j = 1:length(simFineTest_powers)
            for k = countVec
                if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2)
                    line(permute(simFineTest_freqRef(i,j,countVec), [3 2 1])*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9, log10(permute(simFineTest_powerRef(i,j,k), [3 2 1]))*[1 1],...
                        'color', fine_powerCols(j,:), 'linewidth', 2)
                else 
                    if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100
                        plot3(permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, log10(permute(simFineTest_powerRef(i,j,k), [3 2 1])),...
                            '.', 'color', fine_powerCols(j,:), 'markersize', 8)
                    elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) > inactNoiseThresh
                        plot3(permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, log10(permute(simFineTest_powerRef(i,j,k), [3 2 1])),...
                            'o', 'color', fine_powerCols(j,:), 'markersize', 4)
                    end
                end
            end
        end
    end

    ylim([80 120])
    xlim([7 10])

%% Phase 1.4 - test time

if length(simTimeTest_freqs) ~= 1
    warning('1st time test configured for single freq')
end

simTimeTest_inact = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nReps);
simTimeTest_inactRef = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nReps);

simTimeTest_freqRef = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nReps);
simTimeTest_powerRef = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nReps);
simTimeTest_timeRef = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nReps);

% Get matched roots across freq
if length(powerThreshold) == 1
    powerThresholdInterp_powerTest = powerThreshold*ones(length(simTimeTest_freqs),1);
else
    powerThresholdInterp_powerTest = interp1(powerThresholdFreqs, powerThreshold, simTimeTest_freqs, 'linear',0);
    powerThresholdInterp_powerTest(powerThresholdInterp_powerTest == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        simTimeTest_freqs(powerThresholdInterp_powerTest == 0), 'nearest','extrap');
end

for i = 1

    curveVal = curveMax*exp(-(simTimeTest_freqs(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(simTimeTest_powers)

        if simTimeTest_powers(j) > powerThresholdInterp_powerTest(i) 
                    
            if useWeibullPower
                if log10(simTimeTest_powers(j)) >= powerWeibullThreshold
                    powerVal = 1-exp(-powerWeibullAlpha * ...
                        (log10(simTimeTest_powers(j))-powerWeibullThreshold)^powerWeibullBeta);
                else
                   powerVal = 0; 
                end
            else
                powerVal = powerLinearA*(log10(simTimeTest_powers(j)))+powerLinearB;    
            end
        else
            powerVal = 0;
        end
        
        for k = 1:length(simTimeTest_times)
            powerScale = log10(simTimeTest_powers(j))/timeWeibullPowerCenter;
            powerScale = 1 + (powerScale-1)*timeWeibullPowerScale;
            
            if powerScale < 0.01
                powerScale = 0.01;
            end
            
            if log10(simTimeTest_times(k)*powerScale) >= timeWeibullThreshold                   
               timeVal = 1-exp(-timeWeibullAlpha* ...
                        (log10(simTimeTest_times(k)*powerScale)-timeWeibullThreshold)^timeWeibullBeta); 
            else
                timeVal = 0;
            end
            
            totalVal = powerVal * timeVal * curveVal;
            totalVal(totalVal > 100) = 100;     
            totalVal(totalVal < 0) = 0;

            simTimeTest_inactRef(i,j,k,:) = totalVal;
            simTimeTest_inact(i,j,k,:) = totalVal + absStd*randn(nReps,1);
            
            simTimeTest_freqRef(i,j,k,:) = simTimeTest_freqs(i);
            simTimeTest_powerRef(i,j,k,:) = simTimeTest_powers(j);
            simTimeTest_timeRef(i,j,k,:) = simTimeTest_times(k);
        end
    end
end

simTimeTest_inact(simTimeTest_inact > 100) = 100;

% Plot power curve
figure; 
subplot(4,length(simTimeTest_times),1); hold on

for i = 1
    for j = 1:length(simTimeTest_powers)
        for k = 1:length(simTimeTest_times)
            toPlot = setxor(countVec, 1:nReps);

            plot3(log10(permute(simTimeTest_powerRef(i,j,k,toPlot), [4 3 2 1])), permute(simTimeTest_inact(i,j,k,toPlot), [4 3 2 1]), log10(permute(simTimeTest_timeRef(i,j,k,toPlot), [4 3 2 1])),...
                'o', 'markersize', markerSize,'color', time_timeCols(k,:))

            plot3(log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])), permute(simTimeTest_inact(i,j,k,countVec), [4 3 2 1]), log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])),...
                '*', 'markersize', markerSize, 'color', time_timeCols(k,:));
        end
    end
end
plot3(log10(permute(simTimeTest_powerRef(1,:,end,1), [2 1 3 4])), permute(simTimeTest_inactRef(1,:,end,1), [2 1 3 4]), log10(permute(simTimeTest_timeRef(1,:,end,1), [2 1 3 4])),...
    'm', 'linewidth', 2)

plot3(log10(permute(simTimeTest_powerRef(1,:,1,1), [2 1 3 4])), permute(simTimeTest_inactRef(1,:,1,1), [2 1 3 4]), log10(permute(simTimeTest_timeRef(1,:,1,1), [2 1 3 4])),...
    ':m', 'linewidth', 2)

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Plot time curve
subplot(4,length(simTimeTest_times),length(simTimeTest_times)+1); hold on

for i = 1:length(simTimeTest_freqs)
    for j = 1:length(simTimeTest_powers)
        for k = 1:length(simTimeTest_times)
            toPlot = setxor(countVec, 1:nReps);

            plot3(log10(permute(simTimeTest_timeRef(i,j,k,toPlot), [4 3 2 1])), permute(simTimeTest_inact(i,j,k,toPlot), [4 3 2 1]), log10(permute(simTimeTest_powerRef(i,j,k,toPlot), [4 3 2 1])),...
                'o', 'markersize', markerSize,'color', time_powerCols(j,:));

            plot3(log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])), permute(simTimeTest_inact(i,j,k,countVec), [4 3 2 1]), log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])),...
                '*', 'markersize', markerSize, 'color', time_powerCols(j,:));
        end
    end
end

plot3(log10(permute(simTimeTest_timeRef(1,end,:,1), [3 1 2 4])), permute(simTimeTest_inactRef(1,end,:,1), [3 1 2 4]), log10(permute(simTimeTest_powerRef(1,end,:,1), [3 1 2 4])),...
    'm', 'linewidth', 2)

plot3(log10(permute(simTimeTest_timeRef(1,1,:,1), [3 1 2 4])), permute(simTimeTest_inactRef(1,1,:,1), [3 1 2 4]), log10(permute(simTimeTest_powerRef(1,1,:,1), [3 1 2 4])),...
    ':m', 'linewidth', 2)

xlim([-3 3])


inactRatioBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount, length(influenzaSize));

inactFractionBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount, length(influenzaSize));

activeBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount, length(influenzaSize));

totalInact = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount);


inactInds = cell(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times));

maxInactInd = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount); % if more than 1 equal, will just be one 

maxInactBounds = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), length(simTimeTest_times), nrepsCount, 2);

%%% Make a function for this, will use later
for i = 1
    for j = fliplr(1:length(simTimeTest_powers))
        
        for k = fliplr(1:length(simTimeTest_times))
            
            subplot(4,length(simTimeTest_times),2*length(simTimeTest_times)+k); hold on
            plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-b', 'linewidth', 2)
        
            for l = countVec
                tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

                tempSize_samples = sort(tempSize_samples, 'descend');

                tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

                samplesIntact = ones(testCountNum, 1);

                % Start with reference before we figure out noise... 
                numToInact = round(testCountNum*simTimeTest_inactRef(i,j,k,l)/100);

                if numToInact > 1 
                    % Find nearest in sample
                    [~, minInd] = min(abs(tempSize_freqs-simTimeTest_freqs(i)));

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

                            [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-simTimeTest_freqs(i)));

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

                inactRatioBySize(i,j,k,l,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

                inactRatioBySize(i,j,k,l,zeroRatioInds) = 0;

                inactFractionBySize(i,j,k,l,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

                inactFractionBySize(i,j,k,l,zeroFractionInds) = 0;

                activeBySize(i,j,k,l,:) = activeDistDist/testCountNum*100;

                totalInact(i,j,k,l) = (1-sum(activeDistDist)/testCountNum)*100;

                % get inds and ranges on inactivation spectrum

                % find max
                if length(simTimeTest_powers) == j && length(simTimeTest_times) == k
                    referenceRange = find(inactRatioBySize(i,j,k,l,:) == max(inactRatioBySize(i,j,k,l,:)));

                    referenceRange = union(inactInds{i,j,k}, referenceRange);
                elseif length(simTimeTest_powers) == j 
                    referenceRange = inactInds{i, j,k+1};
                    
                else
                    referenceRange = inactInds{i, j+1,k};
                end

                [maxInact, tempInd] = max(inactRatioBySize(i,j,k,l,referenceRange));

                maxInactInd(i,j,k,l) = referenceRange(tempInd);

                % Note that union is used so that this would grow across multiple count replicats
                tempInactInds = referenceRange(inactRatioBySize(i,j,k,l,referenceRange) == maxInact);

                inactInds{i,j,k} = union(tempInactInds, inactInds{i,j,k});

                maxInactBounds(i,j,k,l,:) = influenzaSize(tempInactInds([1 end]));    

                % plotting spectra - counted of inactivated
                subplot(4,length(simTimeTest_times),2*length(simTimeTest_times)+k); hold on

                plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', time_powerCols(j,:))

                % inactivation spectrum
                subplot(4,length(simTimeTest_times),3*length(simTimeTest_times)+k); hold on

                plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,l,:), [5 4 3 2 1]), 'color', time_powerCols(j,:));
                ylim([-50 100])
            end
        end
    end
end

%%% Factored by power, 1st line
    % inact from counted intact summed
    subplot(4,length(simTimeTest_times),4); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot(permute(totalInact(i,j,k,countVec), [4 3 2 1]), permute(simTimeTest_inact(i,j,k,countVec), [4 3 2 1]), 'o', 'color', time_timeCols(k,:), 'markersize', markerSize);
            end
        end
    end
    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]); 

    % Inact from counted both across range
    subplot(4,length(simTimeTest_times),2); hold on

    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot3(log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])), permute(inactFractionBySize(i,j,k,countVec, maxInactInd(i,j,k,countVec)), [5 4 3 2 1]), log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])),...
                    'x', 'markersize', markerSize, 'color', time_timeCols(k,:))
            end
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Inact from counted intact across range
    subplot(4,length(simTimeTest_times),3); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot3(log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])), permute(inactRatioBySize(i,j,k,countVec,maxInactInd(i,j,k,countVec)), [5 4 3 2 1]), log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])), ...
                    'x', 'markersize', markerSize, 'color', time_timeCols(k,:))
            end
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3])

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simTimeTest_times),5); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                for l = countVec
                    if maxInactBounds(i,j,k,l,1) ~= maxInactBounds(i,j,k,l,2)
                        line(log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,l,:), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1]))*[1 1], ...
                            'color', time_timeCols(k,:), 'linewidth', 2)
                    else
                        if inactRatioBySize(i,j,k,l,maxInactInd(i,j,k,l)) == 100
                            plot3(log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1])), permute(maxInactBounds(i,j,k,l,1), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1])), ...
                                '.', 'markersize', 8, 'color', time_timeCols(k,:))
                        elseif inactRatioBySize(i,j,k,l,maxInactInd(i,j,k,l)) > inactNoiseThresh
                            plot3(log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1])), permute(maxInactBounds(i,j,k,l,1), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1])), ...
                                'o', 'markersize', 4, 'color', time_timeCols(k,:))
                        end
                    end
                end
            end
        end
    end
    ylim([80 120])
    xlim([1 3])
    
%%% Factored by time, second line
% inact from counted intact summed
    subplot(4,length(simTimeTest_times),length(simTimeTest_times)+4); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot(permute(totalInact(i,j,k,countVec), [4 3 2 1]), permute(simTimeTest_inact(i,j,k,countVec), [4 3 2 1]), 'o', 'color', time_powerCols(j,:), 'markersize', markerSize);
            end
        end
    end
    line([0 100], [0 100])
    ylim([0 100]); xlim([0 100]); 

    % Inact from counted both across range
    subplot(4,length(simTimeTest_times),length(simTimeTest_times)+2); hold on

    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot3(log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])), permute(inactFractionBySize(i,j,k,countVec, maxInactInd(i,j,k,countVec)), [5 4 3 2 1]), log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])), ...
                    'x', 'markersize', markerSize, 'color', time_powerCols(j,:))
            end
        end
    end
    xlim([-3 3])

    % Inact from counted intact across range
    subplot(4,length(simTimeTest_times),length(simTimeTest_times)+3); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                plot3(log10(permute(simTimeTest_timeRef(i,j,k,countVec), [4 3 2 1])), permute(inactRatioBySize(i,j,k,countVec,maxInactInd(i,j,k,countVec)), [5 4 3 2 1]), log10(permute(simTimeTest_powerRef(i,j,k,countVec), [4 3 2 1])), ...
                    'x', 'markersize', markerSize, 'color', time_powerCols(j,:))
            end
        end
    end
    xlim([-3 3])

    % Plot inactivation width

    %%% Could use thickness or color to indicate number inactivated
    subplot(4,length(simTimeTest_times),length(simTimeTest_times)+5); hold on
    for i = 1
        for j = 1:length(simTimeTest_powers)
            for k = 1:length(simTimeTest_times)
                for l = countVec
                    if maxInactBounds(i,j,k,l,1) ~= maxInactBounds(i,j,k,l,2)
                        line(log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,l,:), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1]))*[1 1], ...
                            'color', time_powerCols(j,:), 'linewidth', 2) 
                    else
                        if inactRatioBySize(i,j,k,l,maxInactInd(i,j,k,l)) == 100
                            plot3(log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1])), permute(maxInactBounds(i,j,k,l,1), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1])), ...
                                '.', 'markersize', 8, 'color', time_powerCols(j,:)) 
                        elseif inactRatioBySize(i,j,k,l,maxInactInd(i,j,k,l)) > inactNoiseThresh
                            plot3(log10(permute(simTimeTest_timeRef(i,j,k,l), [4 3 2 1])), permute(maxInactBounds(i,j,k,l,1), [5 4 3 2 1])*10^9, log10(permute(simTimeTest_powerRef(i,j,k,l), [4 3 2 1])), ...
                                'o', 'markersize', 4, 'color', time_powerCols(j,:)) 
                        end
                    end
                end
            end
        end
    end
    ylim([80 120])
    xlim([-3 3])
    