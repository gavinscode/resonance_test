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
   error('Problem triming - rerun') 
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

measuredSize_resonances = 1./(dataDiameter2017(:,1)/2/10^9)*resSlope;

figure;
subplot(1,2,1); hold on
plot(influenzaSize, influenzaSize_dist);

plot(dataDiameter2017(:,1)/10^9, dataDiameter2017(:,2));

subplot(1,2,2); hold on
plot(influenzSize_resonances, influenzaSize_dist);

plot(measuredSize_resonances, dataDiameter2017(:,2));

%% Parameters for simulations

absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% phenomenlogical parameters
    % For freq
        curveMax = 100;
        curveCenter = 8.5; 
        curveSpread = 2.3;

    % For power - single threhsold
        powerThreshold = 45; %45
        % For power - linear interpolation between points then constant
%         powerThreshold = [25 65]; 
%         powerThresholdFreqs = [7.5 9.5];

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

inactNoiseThresh = 50;

countVec = 1:nrepsCount;

testCountNum = initialCountNum;

% Phase 1.1 - also want to fit function to get peak 
% 2 GHz seperation seems ok to fit curve to inact, but then have less info on individual inact 
simFreqTest_freqs = 1:1:20;
simFreqTest_power = round(10^2.75);
simFreqTest_time = 1000;

freq_freqCols = jet(length(simFreqTest_freqs)); 

% Phase 1.2 - identify power response (minimum effective) at center freq
powerLogRange = [1:0.25:2.75];
%%% Place for symmetry in freq (or for size dist effected?)
% simPowerTest_freqs = [6 7.5 8.5 9.5 11]; 
simPowerTest_freqs = 1:0.5:15.5;
% simPowerTest_freqs = [5.5 6 7 7.5 8.5 9 10 10.5 11.5]; 
% simPowerTest_freqs = [6 8.5 11];

simPowerTest_powers = round(10.^(powerLogRange))
simPowerTest_time = simFreqTest_time;

power_freqCols = winter(length(simPowerTest_freqs));
power_powerCols = cool(length(simPowerTest_powers));

% ranges for plotting
freqRangeFine = 0:0.05:max(simFreqTest_freqs);
powerRangeFine = round(10.^(0:0.005:3)); % 0:1:1000; 

markerSize = 6;
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
            elseif inactRatioBySize(i,j,maxInactInd(i,j)) >= inactNoiseThresh
                plot(simFreqTest_FreqRef(i,j), maxInactBounds(i,j,1)*10^9,...
                    'o', 'markersize', 4, 'color', freq_freqCols(i,:))
            end
        end
    end
end
ylim([80 120])
xlim([1 20])

freqInactRatioBySize = inactRatioBySize;

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
                        elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) >= inactNoiseThresh
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
                    elseif inactRatioBySize(i,j,k,maxInactInd(i,j,k)) >= inactNoiseThresh
                        plot3(permute(simPowerTest_freqRef(i,j,countVec), [3 2 1]), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), ...
                            'o', 'markersize', 4, 'color', power_powerCols(j,:))
                    end
                end
            end
        end
    end
    ylim([80 120])
    
powerInactRatioBySize = inactRatioBySize;

%% Solver using both freq and power scan

removeUnused = 1;

%%% If freq or size spacing used is irregular, will cause interp problems
    % Irregular power spacing ok? (Regular log anyway)

predictionFreqs_init = unique(sort([simFreqTest_freqs simPowerTest_freqs]));
% predictionFreqs = predictionFreqs_init;
predictionFreqs = predictionFreqs_init(1):0.5:predictionFreqs_init(end);

predictionPowers = simPowerTest_powers; %round(10.^(1:0.125:2.75)); 

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

% Remove low lines from power
if removeUnused
    predictInactThresh = 0.1; 
    
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
powerInd = find(predictionPowers == simFreqTest_power);
for i = 1:length(simFreqTest_freqs)
    freqInd = find(simFreqTest_freqs(i) == predictionFreqs);
    
    if ~isempty(freqInd)
        
        meanInacts = mean(freqInactRatioBySize(i,countVec,:),2);
        
        sizeInds = find(meanInacts >= inactNoiseThresh);
        
        [inds] = sub2ind(size(thresholdArray), freqInd*ones(length(sizeInds),1), sizeInds);
        
        thresholdArray(inds) = powerInd;
    end
end

for i = 1:length(simPowerTest_freqs)
    freqInd = find(simPowerTest_freqs(i) == predictionFreqs);
    
    if ~isempty(freqInd)
        for j = fliplr(1:length(simPowerTest_powers))
            
            meanInacts = mean(powerInactRatioBySize(i,j,countVec,:),2);

            sizeInds = find(meanInacts >= inactNoiseThresh);

            [inds] = sub2ind(size(thresholdArray), freqInd*ones(length(sizeInds),1), sizeInds);

            powerInd = find(predictionPowers == simPowerTest_powers(j));
            
            thresholdArray(inds) = powerInd;
        end
    end
end

% Interpolate between treshold lines if extra inds added 
if ~all(ismember(predictionFreqs, predictionFreqs_init))
    % first go along each column and fill zeros from top to bottom - makes sure end lines get connected
    for i = 1:length(influenzaSize)
       lineInds = find(thresholdArray(:,i));
       
       if ~isempty(lineInds)
           lineInds = lineInds(1):lineInds(end);
           
           thresholdArray(lineInds(thresholdArray(lineInds,i) == 0) ,i) = 5;
       end
    end
    
    % then do interpolant - smooths edges a bit
    [sizeRef, freqRef] = meshgrid(1:length(influenzaSize), 1:length(predictionFreqs));
    
    indsCheck = find(thresholdArray == 0);
    indsRef = find(thresholdArray > 0);
    
    tempInterp = scatteredInterpolant(freqRef(indsRef), sizeRef(indsRef), ones(length(indsRef),1), 'nearest', 'none');
    
    tempCheckVals = tempInterp(freqRef(indsCheck), sizeRef(indsCheck));
    
    thresholdArray(indsCheck(tempCheckVals > 0)) = 5;
end

% do first pass to see which sizes would be removed, redo later
toRemove = zeros(length(predictionSizes),1,'logical');

if 1
    sizesToUse = 1:length(influenzaSize); 

    if removeUnused
        for i = 1:length(predictionSizes)
            if all( thresholdArray(:,i) == 0)
                toRemove(i) = 1;
            end
        end

        sizesToUse(toRemove) = [];
    end
end

sizesRemoved = find(toRemove);
predictionCountNum = sum(influenzaSize_dist(sizesToUse));
ignoredCountNum = sum(influenzaSize_dist(toRemove));

% Get reference inactivation for unused 
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

        inactRatioBySize_reference(i,j,:) = (1 - (activeDistDist(:))./...
            (influenzaSize_dist(:)))*100;
        
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

        sizeInds = find(inactRatioBySize_reference(i,j,:) >= inactNoiseThresh);

        [inds] = sub2ind(size(thresholdArray_reference), i*ones(length(sizeInds),1), sizeInds);

        thresholdArray_reference(inds) = j;
    end
end

%%%
% warning('Using reference')
% thresholdArray = thresholdArray_reference;
% predictedInactivation = keptInact;

% Remove from sizes and freqs, as above but on more things
toRemove = zeros(length(predictionSizes),1,'logical');

sizesToUse = 1:length(influenzaSize); 

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
    totalInact(toRemoveTemp,:) = [];
    lostInact(toRemoveTemp,:) = [];
    keptInact(toRemoveTemp,:) = [];
    
    inactRatioBySize_reference(toRemoveTemp,:,:) = [];
    
    if 1
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

        inactRatioBySize_reference(:,:,toRemove) = [];
    end
end

sizesRemoved = find(toRemove);
predictionCountNum = sum(influenzaSize_dist(sizesToUse));
ignoredCountNum = sum(influenzaSize_dist(toRemove));

if any(influenzaSize_dist(sizesToUse)/predictionCountNum*100 < stepPercent)
    warning('small sizes are present')
end

measArray = thresholdArray;

% Indicate lines that are already measured
for i = 1:length(simPowerTest_freqs)
    freqInd = find(simPowerTest_freqs(i) == predictionFreqs);
    measArray(freqInd, thresholdArray(freqInd,:) > 0) = -1;
end    

% Thresholded error for reference
inactivationThreshold_ref = size(predictedInactivation);

for i = 1:length(predictionFreqs)
    for j = 1:length(predictionPowers)

        inds = find(permute(inactRatioBySize_reference(i,j,:), [3 2 1]) >= inactNoiseThresh);
        
        inactivationThreshold_ref(i,j) = sum(predictionSizes_dist(inds));
    end
end

%% Test various interp schemes.
% get inerpolation coords
measInds = find(measArray == -1);
[measX, measY] = ind2sub(size(thresholdArray), measInds);

allInds = find(thresholdArray);
[allX, allY] = ind2sub(size(thresholdArray), allInds);



%%% Note, these distance measures could be weird if irregular spacing on interpolation map

% Do 1d interpolation based on distance to edge
% Get distance map to border
distanceArray = bwdist(~logical(thresholdArray), 'cityblock');

% Need to average value for each distance
distanceVals = unique(distanceArray(measInds));
avgVals = zeros(length(distanceVals),1);

for i = 1:length(distanceVals)
   tempAvgInds = find(distanceArray(measInds) == distanceVals(i));
   
   avgVals(i) = mean(thresholdArray(measInds(tempAvgInds)));
end

cityBlockInterpArray = zeros(size(thresholdArray));

% Nearest and linear basically equal
cityBlockInterpArray(allInds) = interp1(distanceVals, avgVals, distanceArray(allInds), 'nearest', 'extrap');

% above will overwrite measured inds if powers didn't match, so add back
cityBlockInterpArray(measInds) = thresholdArray(measInds);



% Now implement different distance map and vertical option
distanceArrayFromLeft =zeros(size(thresholdArray));

distanceArrayFromRight =zeros(size(thresholdArray));

for i = 1:length(predictionFreqs)
    inds = find(allX == i);
    
    if ~isempty(inds)
        
        yVals = sort(allY(inds));
        
        % Check line does not extend across whole map
        if yVals(1) ~= 1 | yVals(end) ~= size(thresholdArray,2)
            % Check if first value is not on border
            if yVals(1) > 1
                % just increases left to right
                distanceArrayFromLeft(allInds(inds)) = 1:length(inds);
            end

            % Reverse from other side
            if yVals(end) < size(thresholdArray,2)
                 %decreases from left to right   
                 distanceArrayFromRight(allInds(inds)) = fliplr(1:length(inds));
            end
 
        else
           error('Line extends across whole map') 
        end
    end
end

maxOrig = max([ max(distanceArrayFromLeft(:)), max(distanceArrayFromRight(:))]);

%%% Distance map updated to use step back rather than look up

% Always some missing because not placed if on border make look up conversion
overlapInds = find(distanceArrayFromLeft(allInds) > 0 & distanceArrayFromRight(allInds) > 0);
missingLeft = find(distanceArrayFromLeft(allInds) == 0);

referenceVals = unique(distanceArrayFromRight(allInds(overlapInds)));

for i = 1:length(referenceVals)
    % get average for right pixels in overlap
    tempInds = find(distanceArrayFromRight(allInds(overlapInds)) == referenceVals(i));
    
    useValue = mean(distanceArrayFromLeft(allInds(overlapInds(tempInds))));
   
    % find right pixels in missing and place
    tempInds = find(distanceArrayFromRight(allInds(missingLeft)) == referenceVals(i)); 
    
    distanceArrayFromLeft(allInds(missingLeft(tempInds))) = useValue;
end

missedInds = find(distanceArrayFromLeft(allInds(missingLeft)) == 0);

if ~isempty(missedInds)
   error('need to add interpolant') 
end

figure;
subplot(1,2,1)
imshow(distanceArrayFromLeft/maxOrig)

subplot(1,2,2)
imshow(distanceArrayFromRight/maxOrig)

% Manually selected - one is basically just the inverse/rotation of the other
bestDistanceMap = distanceArrayFromLeft;

dualDistInterpArray = zeros(size(thresholdArray));

% Now do interp
if length(simPowerTest_freqs) > 1
    % nearest best for both
    dualDistInterp = scatteredInterpolant(bestDistanceMap(measInds), measX, thresholdArray(measInds), 'nearest', 'nearest');

    dualDistInterpArray(allInds) = dualDistInterp(bestDistanceMap(allInds), allX);
else
    % triangulation fails x included - nearest best
    dualDistInterpArray(allInds) = interp1(bestDistanceMap(measInds), thresholdArray(measInds), distanceArrayFromLeft(allInds), 'nearest', 'extrap');
end

% above will overwrite measured inds if powers didn't match, so add back
dualDistInterpArray(measInds) = thresholdArray(measInds);



% Try 2d scattered interpolant on original data placement
cartInterpArray = zeros(size(thresholdArray));

% note values are put into coordinates, not neccersary if regular spacing on interpolation map
if length(simPowerTest_freqs) > 1
    % nearest best for both
    cartesianInterp = scatteredInterpolant(measX, measY, thresholdArray(measInds), 'nearest', 'nearest');

    cartInterpArray(allInds) = cartesianInterp(allX, allY);
else
    % linear best
    cartInterpArray(allInds) = interp1(measY, thresholdArray(measInds), allY, 'nearest', 'extrap');
end

% plot threshold maps and diffs
figure; 
subplot(4,2,1); imshow(thresholdArray/length(predictionPowers));
title('Input')

subplot(4,2,2); imshow(thresholdArray_reference/length(predictionPowers));
title('Reference')

subplot(4,2,3); imshow(cityBlockInterpArray/length(predictionPowers));
title('Initial interp')

subplot(4,2,4); imshow(abs(cityBlockInterpArray - thresholdArray_reference)/length(predictionPowers));
title('Dif start interp to 2d interp')

sseCityBlock_Thresh = sum((cityBlockInterpArray(:) - thresholdArray_reference(:)).^2);

subplot(4,2,5); imshow(cartInterpArray/length(predictionPowers));
title('2d interp')

subplot(4,2,6); imshow(abs( cartInterpArray - thresholdArray_reference)/length(predictionPowers));
title('Dif 2d interp to ref')

sse2D_Thresh = sum((cartInterpArray(:) - thresholdArray_reference(:)).^2);

subplot(4,2,7); imshow(dualDistInterpArray/length(predictionPowers));
title('Dual dist interp')

subplot(4,2,8); imshow(abs( dualDistInterpArray - thresholdArray_reference)/length(predictionPowers));
title('Dif dual dist interp to ref')

sseDual_Thresh = sum((dualDistInterpArray(:) - thresholdArray_reference(:)).^2);

errorThresh = [sseCityBlock_Thresh, sse2D_Thresh, sseDual_Thresh]

% Plot error on inactivaiton
tempPoints = find(thresholdArray_reference);

% firstly against full inact
fun = @(x)inactivationError(x, (thresholdArray), predictionSizes_dist, ... 
    predictedInactivation, tempPoints, 1:length(predictionPowers)); 

[sseCityBlock_inact, ~, cityBlockInactMap] = fun(cityBlockInterpArray(tempPoints));

[sseCart_inact, ~, cartInactMap] = fun(cartInterpArray(tempPoints));

[sseDual_inact, ~, dualInactMap] = fun(dualDistInterpArray(tempPoints));

errorFull = [sseCityBlock_inact, sseCart_inact, sseDual_inact]

figure;
subplot(4,4,2)
imshow(predictedInactivation/100)
title('Ref full inact')

subplot(4,4,3)
imshow(inactivationThreshold_ref/100)
title('Ref threshold inact')

subplot(4,4,5)
imshow(cityBlockInactMap/100)
title('City block inact')

subplot(4,4,6)
imshow(abs(predictedInactivation-cityBlockInactMap)/10)
title('Dif full to city')

subplot(4,4,9)
imshow(cartInactMap/100)
title('2D  inact')

subplot(4,4,10)
imshow(abs(predictedInactivation-cartInactMap)/10)
title('Dif full to 2d')

subplot(4,4,13)
imshow(dualInactMap/100)
title('Dual inact')

subplot(4,4,14)
imshow(abs(predictedInactivation-dualInactMap)/10)
title('Dif full to dual')

subplot(4,4,8); hold on
plot(predictionFreqs, sum(abs(predictedInactivation-cityBlockInactMap),2),'g')
plot(predictionFreqs, sum(abs(predictedInactivation-cartInactMap),2),'r')
plot(predictionFreqs, sum(abs(predictedInactivation-dualInactMap),2),'b')

% second on threshold (overrights previous)
fun = @(x)inactivationError(x, (thresholdArray), predictionSizes_dist, ... 
    inactivationThreshold_ref, tempPoints, 1:length(predictionPowers)); 

[sseCityBlock_inact, ~, cityBlockInactMap] = fun(cityBlockInterpArray(tempPoints));

[sseCart_inact, ~, cartInactMap] = fun(cartInterpArray(tempPoints));

[sseDual_inact, ~, dualInactMap] = fun(dualDistInterpArray(tempPoints));

errorFullThresh = [sseCityBlock_inact, sseCart_inact, sseDual_inact]

subplot(4,4,7)
imshow(abs(inactivationThreshold_ref-cityBlockInactMap)/10)
title('Dif thresh to city block')

subplot(4,4,11)
imshow(abs(inactivationThreshold_ref-cartInactMap)/10)
title('Dif thresh to cart')

subplot(4,4,15)
imshow(abs(inactivationThreshold_ref-dualInactMap)/10)
title('Dif thresh to dual')

subplot(4,4,16); hold on
plot(predictionFreqs, sum(abs(inactivationThreshold_ref-cityBlockInactMap),2),'g')
plot(predictionFreqs, sum(abs(inactivationThreshold_ref-cartInactMap),2),'r')
plot(predictionFreqs, sum(abs(inactivationThreshold_ref-dualInactMap),2),'b')

%% Test 3D interp - fml

%%% Take as given from previous map, may wish to refine boundaries later

% Test build map using full base and centre line...
[sizeGrid, freqGrid, powerGrid] = meshgrid(1:length(sizesToUse), 1:length(predictionFreqs), 1:length(predictionPowers));

referenceVol = permute(inactRatioBySize_reference, [1 3 2]);

inactivationData = referenceVol;

testThresh = 10;

% Clip low values for now
% referenceVol(referenceVol < inactNoiseThresh) = 0;

% build exisiting arrays into volumes
measVol = zeros(size(referenceVol));
distVol = zeros(size(referenceVol));

includedVol = referenceVol;
includedVol(includedVol > 0) = 1;

for i = fliplr(1:length(predictionPowers))
    
    tempArray = thresholdArray*0;
    
    if i == length(predictionPowers)
        compFreqs = predictionFreqs_init;
    else
        compFreqs = simPowerTest_freqs;
    end
    
    for j = 1:length(compFreqs)
        freqInd = find(compFreqs(j) == predictionFreqs);
        tempArray(freqInd, referenceVol(freqInd,:,i) > testThresh) = -1;
    end    

    measVol(:,:,i) = tempArray;
    
    seedArray =  measVol(:,:,i);
    % If first, do interpolation across
    if i == length(predictionPowers)
        if ~all(ismember(predictionFreqs, predictionFreqs_init))
            % first go along each column and fill zeros from top to bottom - makes sure end lines get connected
            for j = 1:length(sizesToUse)
               lineInds = find(seedArray(:,j));

               if ~isempty(lineInds)
                   lineInds = lineInds(1):lineInds(end);

                   seedArray(lineInds(seedArray(lineInds,j) == 0) ,j) = -1;
               end
            end

            % then do interpolant - smooths edges a bit
            [sRef, fRef] = meshgrid(1:length(sizesToUse), 1:length(predictionFreqs));

            indsCheck = find(seedArray == 0);
            indsRef = find(seedArray > 0);

            tempInterp = scatteredInterpolant(fRef(indsRef), sRef(indsRef), ones(length(indsRef),1), 'nearest', 'none');

            tempCheckVals = tempInterp(fRef(indsCheck), sRef(indsCheck));

            seedArray(indsCheck(tempCheckVals > 0)) = -1;   
        end
        
        baseSeedArray = seedArray;
    end
    
    % Do dist map as for 2D
    tempArrayFromLeft = zeros(size(thresholdArray));
    tempArrayFromRight = zeros(size(thresholdArray));

    levelInds = find(seedArray == -1);
    [levelX, levelY] = ind2sub(size(thresholdArray), levelInds);
    
    maxOpenLeft = [];
    maxOpenRight = [];
    
    for j = 1:length(predictionFreqs)
        inds = find(levelX == j);

        if ~isempty(inds)

            yVals = sort(levelY(inds));

            % Check line does not extend across whole map
            if yVals(1) ~= 1 | yVals(end) ~= size(thresholdArray,2)
                % Check if first value is not on border
                if yVals(1) > 1
                    % just increases left to right
                    tempArrayFromLeft(levelInds(inds)) = 1:length(inds);
                    
                    if yVals(end) < size(thresholdArray,2)
                        maxOpenLeft = [maxOpenLeft length(inds)];
                    end
                end

                % Reverse from other side
                %%% Not really needed now
                if yVals(end) < size(thresholdArray,2)
                     %decreases from left to right   
                     tempArrayFromRight(levelInds(inds)) = fliplr(1:length(inds));
                     
                     if yVals(end) > 1
                        maxOpenRight = [maxOpenRight length(inds)];
                     end
                end

            else
               error('Line extends across whole map') 
            end
        end
    end

    missingLeft = find(tempArrayFromLeft(levelInds) == 0);
    
    % Always some missing because not placed if on border make look up conversion
%     overlapInds = find(tempArrayFromLeft(levelInds) > 0 & tempArrayFromRight(levelInds) > 0);

%     referenceVals = unique(tempArrayFromRight(levelInds(overlapInds)));   
%     
%     for j = 1:length(referenceVals)
%         % get average for right pixels in overlap
%         tempInds = find(tempArrayFromRight(levelInds(overlapInds)) == referenceVals(j));
% 
%         useValue = mean(tempArrayFromLeft(levelInds(overlapInds(tempInds))));
% 
%         % find right pixels in missing and place
%         tempInds = find(tempArrayFromRight(levelInds(missingLeft)) == referenceVals(j));
% 
%         tempArrayFromLeft(levelInds(missingLeft(tempInds))) = useValue;
%     end

    if all(maxOpenLeft ~= maxOpenLeft(1))
       error('Need to revise') 
    end

    maxOpenLeft = maxOpenLeft(1);
    
    % Changed to using step back instead of look up table
    for j = 1:length(predictionFreqs)
        inds = find(levelX(missingLeft) == j);
        
        if ~isempty(inds)
            yVals = sort(levelY(missingLeft(inds)));

            % Check line does not extend across whole map
            if yVals(1) ~= 1 | yVals(end) ~= size(thresholdArray,2)
                
                tempArrayFromLeft(levelInds(missingLeft(inds))) =  (maxOpenLeft-length(inds)+1):maxOpenLeft;
                
                % if border is low, step it up
                if tempArrayFromLeft(levelInds(missingLeft(inds(1)))) <= 0
                    tempArrayFromLeft(levelInds(missingLeft(inds))) = tempArrayFromLeft(levelInds(missingLeft(inds))) + 1 - tempArrayFromLeft(levelInds(missingLeft(inds(1))));
                end
            end
        end
    end
    
    % On border - may need to fix in future
    missedInds = find(tempArrayFromLeft(levelInds(missingLeft)) == 0)
    
    % If not base, interpolate up as well
    if 0; %length(predictionPowers) ~= i
        baseInds = find(baseSeedArray == -1);
        baseArrayFromLeft = distVol(:,:,end);
        
        overlapInds = find(tempArrayFromLeft(baseInds) > 0 & baseArrayFromLeft(baseInds) > 0);
        missingLeft = find(tempArrayFromLeft(baseInds) == 0);
        
        referenceVals = unique(baseArrayFromLeft(baseInds(overlapInds)));   
    
        for j = 1:length(referenceVals)
            % get average for right pixels in overlap
%             tempInds = find(baseArrayFromLeft(baseInds(overlapInds)) == referenceVals(j));
% 
%             useValue = mean(tempArrayFromLeft(baseInds(overlapInds(tempInds))));

            % find right pixels in missing and place
            tempInds = find(baseArrayFromLeft(baseInds(missingLeft)) == referenceVals(j));

            tempArrayFromLeft(baseInds(missingLeft(tempInds))) = 1;
        end
        
    end
    
    distVol(:,:,i) = tempArrayFromLeft; %distanceArrayFromLeft;
end

% To test
% measVol = referenceVol;
% measVol(measVol > 0) = -1;

% measVol(measVol == -1 & distVol == 0) = 0;

% Dist vol still copied in fully...
measInds = find(measVol == -1);

allInds = find(distVol);
% allInds = find(referenceVol > 0);

find(distVol == 0 & referenceVol > testThresh)

interpVol = zeros(size(referenceVol));

% Now do interp
if length(simPowerTest_freqs) > 1
    % nearest best for both
    distInterp = scatteredInterpolant(distVol(measInds), freqGrid(measInds), powerGrid(measInds), referenceVol(measInds), 'linear', 'none');

    interpVol(allInds) = distInterp(distVol(allInds), freqGrid(allInds), powerGrid(allInds));
else
    % triangulation fails x included - nearest best
    distInterp = scatteredInterpolant(distVol(measInds), powerGrid(measInds), referenceVol(measInds), 'nearest', 'nearest');

    interpVol(allInds) = distInterp(distVol(allInds), powerGrid(allInds));
end

% above will overwrite measured inds if powers didn't match, so add back
interpVol(measInds) = referenceVol(measInds);

interpVol(isnan(interpVol(:))) = 0;

figure;
for i = 1:5
    subplot(5,3,i*3-2)
    imshow(referenceVol(:,:,i)/100)
    
    subplot(5,3,i*3-1)
    imshow(interpVol(:,:,i)/100)
    
    subplot(5,3,i*3)
    imshow(distVol(:,:,i))
end

% get inactivation
interpInactArray = size(predictedInactivation);

for i = 1:length(predictionFreqs)
    for j = 1:length(predictionPowers)

        interpInactArray(i,j) = sum(predictionSizes_dist/100.*interpVol(i,:,j));
    end
end

% plot error
figure;
subplot(4,4,2)
imshow(predictedInactivation/100)
title('Ref full inact')

subplot(4,4,3)
imshow(keptInact*predictionCountNum/initialCountNum/100)
title('Ref kept inact')

subplot(4,4,5)
imshow(interpInactArray/100)
title('From volume')

subplot(4,4,6)
imshow(abs(predictedInactivation-interpInactArray)/10)
title('Dif full to vol')

subplot(4,4,7)
imshow(abs(keptInact*predictionCountNum/initialCountNum-interpInactArray)/10)
title('Dif kept to vol')

subplot(4,4,8); hold on
plot(predictionFreqs, sum(abs(predictedInactivation-interpInactArray),2),'g')

subplot(4,4,16); hold on
plot(predictionFreqs, sum(abs(keptInact*predictionCountNum/initialCountNum-interpInactArray),2),'g')

sseVol_full = sum((predictedInactivation(:)-interpInactArray(:)).^2);
sseVol_thresh = sum((keptInact(:)*predictionCountNum/initialCountNum-interpInactArray(:)).^2);

[sseVol_full, sseVol_thresh]