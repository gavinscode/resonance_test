clear; close all; clc;
%% Set parameters for simulations

% Distribution parameters
    influenzaSize_mean = 100*10^-9;
    influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

    initialCountNum = 5000; 
    
    %%% May be some optimum way to set size to get best SNR later.
        % How reliably can size be determined from TEM? surely not < 1 nm... 
    sizeStep = 2; % in nm
    sizeRange = [70 130];


    % lowering this reduces error on ref, but will give noise on real
    sizePercent_thresh = 2; % Note that nothing below this percent will get included in size maps
    
    %%% In making final plots this thresh probal interacts with predict and inact thresholds
    
    %%% This should be determined by matching peak freq to size, change later
    influenzaVl = 1486; % Just matched in Saviot tool
    vLtoVtRatio = 0.5;
    influenzaVt = influenzaVl*vLtoVtRatio;

% phenomenlogical parameters
    % For freq
        curveMax = 100;
        curveCenter = 8.5; 
        curveSpread = 2.3;

    % For power - single threhsold
        %%% based on frequency, should adjust to be size dependent..
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
        
% Test parameters
    % Inactivation parameters
        absStd= 5; % std on inactivation, absolute, not relative values (lower SNR on low inactiviation)
        
        nReps = 1; % Difference is legacy from including plaque assay
        nrepsCount = nReps; % for EM - will be match to one plaque assay, but independent in simulation
             countVec = 1:nrepsCount;
        
        inactNoiseThresh = 50; % on EM

        testCountNum = -1; % if -1, will set equal to reference after distribution clipped

% Frequency scan
    simFreqTest_freqs = 1:1:20;
    freqPowerLog = 2.8;
    simFreqTest_power = round(10^freqPowerLog); % set lower than expected saturation

    freq_freqCols = jet(length(simFreqTest_freqs)); 

% Power scanning

    % inital power range set - reduce range and finer steps 
    powerLog = [1:0.2:freqPowerLog];
    simPowerTest_powers = round(10.^(powerLog))

    %%% currently freqs set manually for testing, but should set frequency to center in 1st pass and then minimize error on later passes...
%     simPowerTest_freqs = 8.5;
        
    % all freqs
%     simPowerTest_freqs = 1:0.5:15.5;

    % common progression to minimize error
  simPowerTest_freqs = [6.5 8.5 10.5];
%   simPowerTest_freqs = [6.5 7.5 8.5 10.5 11.5];
    
    % problem combos
%   simPowerTest_freqs = [4.5 8.5 12.5]; % wide
%   simPowerTest_freqs =  [8.5 11.5]; % [5.5 8.5]; % pair
%   simPowerTest_freqs = [5.5 8.5 12.5]; % [4.5 8.5 11.5]; % assym

    power_freqCols = winter(length(simPowerTest_freqs));
    power_powerCols = cool(length(simPowerTest_powers));

% For reference calculation
    referenceFreqRange = min(simFreqTest_freqs):0.5:max(simFreqTest_freqs);
    referencePowerRange = round(10.^(min(powerLog):0.2:max(powerLog)));
        
    % Freq and power divisions used in test should overlap to reference
    if ~all(ismember(simFreqTest_freqs, referenceFreqRange))
        error('Test freqs missing from reference calc')
    end
    
    if ~all(ismember(simPowerTest_powers, referencePowerRange))
        error('Test powers missing from reference calc')
    end
    
% For prediction
    predFreqStep = 0.5;
    
    predictInactThresh = 0.5; % ignore if total inact is under this percentage 
    
% For interpolation    

    % How high can this be taken without changing results?
        %%% 10 v 20 changes increases when using all rows
    % 0 for ideal data, but will cause error if missing rows...
    interpInactThresh = 10; % don't include if indvidual size inact is under this percentage
    
    % distInterp = 'nearest'; % best if wide spacing...
    distInterpType = 'linear'; % best for most

% For plotting    
    % ranges for plotting fitted curves
    freqRangeFine = 0:0.05:max(simFreqTest_freqs);
    powerRangeFine = round(10.^(0:0.05:3)); % 0:1:1000; 

    markerSize = 6;    
    
%% Load digitized datasets
compare_influenza_data
close all
%% generate parameters for size - would be from SEM

influenzaSize_samples = randn(initialCountNum,1)*influenzaSize_std + influenzaSize_mean;

[influenzaSize_dist, influenzaSize] = hist(influenzaSize_samples,(sizeRange(1):sizeStep:sizeRange(2))/10^9);

%%% Skipping for now
% Trim distribution and samples
% tempInd = find(influenzaSize_dist > sum(influenzaSize_dist)*distPercent_thresh/100);
% 
% %%% Note that bin centers are specified
% influenzaSize_samples(influenzaSize_samples < mean(influenzaSize(tempInd(1)-1:tempInd(1))) | ...
%     influenzaSize_samples > mean(influenzaSize(tempInd(end):tempInd(end)+1))) = [];
% 
% influenzaSize(influenzaSize_dist < sum(influenzaSize_dist)*distPercent_thresh/100) = [];
% influenzaSize_dist(influenzaSize_dist < sum(influenzaSize_dist)*distPercent_thresh/100) = [];
% 
% initialCountNum = sum(influenzaSize_dist);

if testCountNum == -1
    testCountNum = initialCountNum;
end

if length(influenzaSize_samples) ~= initialCountNum
   error('Problem triming - rerun') 
end

[~, maxSizeInd] = max(influenzaSize_dist);

sizes2Use = find(influenzaSize_dist/initialCountNum*100 >= sizePercent_thresh);
sizes2Remove = find(influenzaSize_dist/initialCountNum*100 < sizePercent_thresh);
sizes2RemoveLow = sizes2Remove(sizes2Remove < sizes2Use(1));
sizes2RemoveHigh = sizes2Remove(sizes2Remove > sizes2Use(end));
%% get resonance to use for simulation - not known in experiment 
tempResSmall = calcualtesphereresonance(influenzaSize(1)/2, ...
                'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
            
tempResLarge = calcualtesphereresonance(influenzaSize(end)/2, ...
    'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;

% Note resonance of small is higher than large - get slope
resSlope_sim = (tempResSmall-tempResLarge)/(1/(influenzaSize(1)/2)-1/(influenzaSize(end)/2));

%% Calculate reference inactivation for phenom and size distribution

% Get matched power across freq
if length(powerThreshold) == 1
    powerThresholdInterp_ref = powerThreshold*ones(length(referenceFreqRange),1);
else
    % Interp if more than one
    powerThresholdInterp_ref = interp1(powerThresholdFreqs, powerThreshold, referenceFreqRange, 'linear',0);
    powerThresholdInterp_ref(powerThresholdInterp_ref == 0) = interp1(powerThresholdFreqs, powerThreshold, ...
        referenceFreqRange(powerThresholdInterp_ref == 0), 'nearest','extrap');
end

[fullPowerRef, fullFreqRef] = meshgrid(referencePowerRange, referenceFreqRange);

referenceInact = zeros(length(referenceFreqRange),length(referencePowerRange));

for i = 1:length(referenceFreqRange)
    curveVal = curveMax*exp(-(referenceFreqRange(i)-curveCenter).^2/(2*curveSpread^2));
    
    for j = 1:length(referencePowerRange)
        if referencePowerRange(j) > powerThresholdInterp_ref(i)
            if useWeibullPower
                if log10(referencePowerRange(j)) >= powerWeibullThreshold
                    powerVal = 1-exp(-powerWeibullAlpha * ...
                        (log10(referencePowerRange(j))-powerWeibullThreshold)^powerWeibullBeta);
                else
                   powerVal = 0; 
                end
            else
                powerVal = powerLinearA*(log10(referencePowerRange(j)))+powerLinearB;    
            end
        else
            powerVal = 0;
        end

        totalVal = powerVal * curveVal;
        totalVal(totalVal > 100) = 100;     
        totalVal(totalVal < 0) = 0;

        referenceInact(i,j) = totalVal;
    end
end

% Now get across distribution
referenceInactRatio = zeros(length(referenceFreqRange), length(referencePowerRange), length(influenzaSize));

for i = 1:length(referenceFreqRange)
    for j = fliplr(1:length(referencePowerRange))
        tempSize_samples = influenzaSize_samples;

        tempSize_samples = sort(tempSize_samples, 'descend');

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope_sim;

        samplesIntact = ones(initialCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(initialCountNum*referenceInact(i,j)/100);

        if numToInact > 1 
            % Find nearest in sample
            [~, minInd] = min(abs(tempSize_freqs-referenceFreqRange(i)));

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

                    [~, tempInd] = min(abs(tempSize_freqs(samplesIntact == 1)-referenceFreqRange(i)));

                    minInd = indRef(tempInd);
                end

                samplesIntact(minInd) = 0;

                numInactivated = numInactivated + 1;
            end
        end

        % Get hist of both
        activeDistDist = hist(tempSize_samples(samplesIntact == 1), influenzaSize);

        referenceInactRatio(i,j,:) = (1 - (activeDistDist(:))./...
            (influenzaSize_dist(:)))*100;
    end
end

referenceInactRatio(isnan(referenceInactRatio)) = 0;
%% Frequency scan

freqScanIndexToRef = zeros(length(simFreqTest_freqs), 1);

for i = 1:length(simFreqTest_freqs)
    freqScanIndexToRef(i) = find(fullFreqRef == simFreqTest_freqs(i) & fullPowerRef == simFreqTest_power);
end

%%% Need to save this for later and write into combined array...

inactRatioToOrigBySize = zeros(length(simFreqTest_freqs), length(influenzaSize), nrepsCount); % From ratio to original dist
inactRatioEachBySize = zeros(length(simFreqTest_freqs), length(influenzaSize), nrepsCount); % From ratio on each (if inact countable)

simFreqTest_inact_EM = zeros(length(simFreqTest_freqs), nrepsCount);

maxInactInd = zeros(length(simFreqTest_freqs), nrepsCount); % if more than 1 equal, will just be one 

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

        tempSize_freqs = 1./(tempSize_samples/2)*resSlope_sim;

        samplesIntact = ones(testCountNum, 1);

        % Start with reference before we figure out noise... 
        numToInact = round(testCountNum*referenceInact(freqScanIndexToRef(i))/100);
        
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
        inactiveDist = hist(tempSize_samples(samplesIntact == 0), influenzaSize);

        activeDistDist = hist(tempSize_samples(samplesIntact == 1), influenzaSize);

        % Strictly speaking, to match equivelent of usual Ct/C0 doesn't require distribution of inactivated
        %%% Should know proper ratio, which may require fitting sides of active distribution
        inactRatioToOrigBySize(i,:,j) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

        inactRatioEachBySize(i,:,j) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

        simFreqTest_inact_EM(i,j) = (1-sum(activeDistDist)/testCountNum)*100;
        
        tempInactRatio = inactRatioToOrigBySize(i,sizes2Use,j);
        
        % get inds and ranges on inactivation spectrum - find max
        referenceRange = find(tempInactRatio == max(tempInactRatio));

        [maxInact, tempInd] = max(tempInactRatio(referenceRange));

        maxInactInd(i,j) = sizes2Use(referenceRange(tempInd)); 


        % plotting spectra - counted of inactivated
        subplot(3,5,5+freqPlot); hold on

        plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', freq_freqCols(i,:))

        % inactivation spectrum
        subplot(3,5,2*5+freqPlot); hold on

        plot(influenzaSize(sizes2Use)*10^9, inactRatioToOrigBySize(i,sizes2Use,j), 'color', freq_freqCols(i,:));
        
        plot(influenzaSize(sizes2RemoveLow)*10^9, inactRatioToOrigBySize(i,sizes2RemoveLow,j), ':', 'color', freq_freqCols(i,:));
        plot(influenzaSize(sizes2RemoveHigh)*10^9, inactRatioToOrigBySize(i,sizes2RemoveHigh,j), ':', 'color', freq_freqCols(i,:));
        ylim([-10 110])
    end
end

% fit curve to inactivation
opts = fitoptions('gauss1', 'Lower', [0 1 0], 'Upper', [100 20 10]);

freqCurve = fit(simFreqTest_freqs(:), simFreqTest_inact_EM(:), 'gauss1', opts);

freq_curveInact = freqCurve(freqRangeFine);
freq_coefBound = confint(freqCurve);

freq_curveConfidence = predint(freqCurve, freqRangeFine, 0.95, 'Functional');
freq_curveObserved = predint(freqCurve, freqRangeFine, 0.95, 'obs');

%%% Set flag for when to use this
%simPowerTest_freqs = freqCurve.b1;

% general plotting
subplot(3,5,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
for i = 1:length(simFreqTest_freqs)
    plot(simFreqTest_freqs(i), simFreqTest_inact_EM(i), ...
        'x', 'markersize', markerSize, 'color', freq_freqCols(i,:));
end
xlim([1 20]); ylim([-10 110])

plot(freqRangeFine, freq_curveInact, 'r-', 'linewidth', 2);
plot(freqRangeFine, freq_curveConfidence, 'r--', 'linewidth', 2);

% Inact from counted both across range
subplot(3,5,2); hold on

for i = 1:length(simFreqTest_freqs)
    plot(simFreqTest_freqs(i), inactRatioEachBySize(i,maxInactInd(i,countVec), countVec),...
        'x', 'markersize', markerSize, 'color', freq_freqCols(i,:))
end
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
xlim([1 20]);  ylim([-10 110])

% Inact from counted intact across range
subplot(3,5,3); hold on
for i = 1:length(simFreqTest_freqs)
        plot(simFreqTest_freqs(i), inactRatioToOrigBySize(i,maxInactInd(i,countVec),countVec),...
            'x', 'markersize', markerSize, 'color', freq_freqCols(i,:))
end
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
xlim([1 20]);  ylim([-10 110])

%% Calculate frequency dist of sizes

%%% Using preset VL, change to fit VL between measured peak freq and fitted peak size
    %%% Can refine this using EM to get precise size inactivated at given freq in power scan

% Get resonance for largest size
tempResSmall = calcualtesphereresonance(influenzaSize(1)/2, ...
                'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;
            
tempResLarge = calcualtesphereresonance(influenzaSize(end)/2, ...
    'sph', 1, 0, influenzaVl, influenzaVt, 10^9, 10^6, 0)/10^9;

% Note resonance of small is higher than large - get slope
resSlope = (tempResSmall-tempResLarge)/(1/(influenzaSize(1)/2)-1/(influenzaSize(end)/2));

% Inverse sizes
influenzSize_resonances = 1./(influenzaSize/2)*resSlope;
measuredSize_resonances = 1./(dataDiameter2017(:,1)/2/10^9)*resSlope;

figure;
subplot(1,2,1); hold on
plot(influenzaSize, influenzaSize_dist/sum(influenzaSize_dist));
plot(dataDiameter2017(:,1)/10^9, dataDiameter2017(:,2)/sum(dataDiameter2017(:,2)));
title('Simulated vs measured size dist')

subplot(1,2,2); hold on
plot(influenzSize_resonances, influenzaSize_dist/sum(influenzaSize_dist));
plot(measuredSize_resonances, dataDiameter2017(:,2)/sum(dataDiameter2017(:,2)));
plot(dataAbs2016(:,1), dataAbs2016(:,2)/100, 'k', 'linewidth',2)

title('Simulated vs measured freq dist')

%% Phase 1.2 - power

powerScanIndexToRef = zeros(length(simPowerTest_freqs), length(simPowerTest_powers));

if ~all(ismember(simPowerTest_freqs, referenceFreqRange))
    error('Power freq is missing from reference calc')
end

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        powerScanIndexToRef(i,j) = find(fullFreqRef == simPowerTest_freqs(i) & fullPowerRef == simPowerTest_powers(j));
    end
end

inactRatioToOrigBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), length(influenzaSize), nrepsCount);
inactRatioEachBySize = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), length(influenzaSize), nrepsCount);

simPowerTest_inact_EM = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount);

inactInds = cell(length(simPowerTest_freqs), length(simPowerTest_powers));
maxInactInd = zeros(length(simPowerTest_freqs), length(simPowerTest_powers), nrepsCount); % if more than 1 equal, will just be one 

figure; 

%%% Make a function for this, will use later
for i = 1:length(simPowerTest_freqs)
    for j = fliplr(1:length(simPowerTest_powers))

        subplot(4,length(simPowerTest_powers),2*length(simPowerTest_powers)+j); hold on
        plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-k', 'linewidth', 2)

        for k = countVec
            tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

            tempSize_samples = sort(tempSize_samples, 'descend');

            tempSize_freqs = 1./(tempSize_samples/2)*resSlope_sim;

            samplesIntact = ones(testCountNum, 1);

            % Start with reference before we figure out noise... 
            numToInact = round(testCountNum*referenceInact(powerScanIndexToRef(i,j))/100);

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

            
            inactRatioToOrigBySize(i,j,:,k) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

            inactRatioEachBySize(i,j,:,k) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

            simPowerTest_inact_EM(i,j,k) = (1-sum(activeDistDist)/testCountNum)*100;
            
            tempInactRatio = inactRatioToOrigBySize(i,j,sizes2Use,k);
            
            % get inds and ranges on inactivation spectrum
            % find max
            if length(simPowerTest_powers) == j
                referenceRange = find(tempInactRatio == max(tempInactRatio));
                
                referenceRange = union(inactInds{i,j}, referenceRange);
            else
                referenceRange = inactInds{i, j+1};
            end
            
            [maxInact, tempInd] = max(tempInactRatio(referenceRange));
            
            maxInactInd(i,j,k) = sizes2Use(referenceRange(tempInd));
            
            % Note that union is used so that this would grow across multiple count replicates
            tempInactInds = referenceRange(tempInactRatio(referenceRange) == maxInact);
            
            inactInds{i,j} = union(tempInactInds, inactInds{i,j});
            
            % plotting spectra - counted of inactivated
            subplot(4,length(simPowerTest_powers),2*length(simPowerTest_powers)+j); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', power_freqCols(i,:))

            % inactivation spectrum
            subplot(4,length(simPowerTest_powers),3*length(simPowerTest_powers)+j); hold on

            plot(influenzaSize(sizes2Use)*10^9, permute(inactRatioToOrigBySize(i,j,sizes2Use,k), [3 2 1 4]), 'color', power_freqCols(i,:))
            plot(influenzaSize(sizes2RemoveLow)*10^9, permute(inactRatioToOrigBySize(i,j,sizes2RemoveLow,k), [3 2 1 4]), ':', 'color', power_freqCols(i,:))
            plot(influenzaSize(sizes2RemoveHigh)*10^9, permute(inactRatioToOrigBySize(i,j,sizes2RemoveHigh,k), [3 2 1 4]), ':', 'color', power_freqCols(i,:))
            ylim([-10 110])
        end
    end
end

%%% Need to fit phenomological power eqn


% General plotting
% Plot power curve
subplot(4,length(simPowerTest_powers),1); hold on

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)

        plot3(log10(simPowerTest_powers(j)), permute(simPowerTest_inact_EM(i,j,countVec), [3 2 1]), simPowerTest_freqs(i), ...
            'x', 'markersize', markerSize, 'color', power_freqCols(i,:))
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3]); ylim([-10 110])

% Plot freq curve
subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+1); hold on

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        plot3(simPowerTest_freqs(i), permute(simPowerTest_inact_EM(i,j,countVec), [3 2 1]), log10(simPowerTest_powers(j)), ...
            'x', 'markersize', markerSize, 'color', power_powerCols(j,:))
    end
end

xlim([5 12]); ylim([-10 110])
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)

%%% Plot factored by power
    % Inact from counted both across range
    subplot(4,length(simPowerTest_powers),2); hold on

    for i = 1:length(simPowerTest_freqs) 
        for j = 1:length(simPowerTest_powers)
            plot3(log10(simPowerTest_powers(j)), permute(inactRatioEachBySize(i,j,maxInactInd(i,j,countVec),countVec), [4 3 2 1]), simPowerTest_freqs(i), ...
                'x', 'markersize', markerSize, 'color', power_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3]); ylim([-10 110])

    % Inact from counted intact across range
    subplot(4,length(simPowerTest_powers),3); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers)
            plot3(log10(simPowerTest_powers(j)), permute(inactRatioToOrigBySize(i,j,maxInactInd(i,j,countVec),countVec), [4 3 2 1]), simPowerTest_freqs(i), ...
                'x', 'markersize', markerSize, 'color', power_freqCols(i,:))
        end
    end
    plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
    xlim([1 3]);  ylim([-10 110])
    
%%% Plot factored by freq
    % Inact from counted both across range
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+2); hold on

    for i = 1:length(simPowerTest_freqs) 
        for j = 1:length(simPowerTest_powers) 
            plot3(simPowerTest_freqs(i), permute(inactRatioEachBySize(i,j, maxInactInd(i,j,countVec),countVec), [4 3 2 1]), log10(simPowerTest_powers(j)), ...
                'x', 'markersize', markerSize, 'color', power_powerCols(j,:))
        end
    end
    xlim([5 12]);  ylim([-10 110])
    plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
    
    % Inact from counted intact across range
    subplot(4,length(simPowerTest_powers),length(simPowerTest_powers)+3); hold on
    for i = 1:length(simPowerTest_freqs)
        for j = 1:length(simPowerTest_powers) 
            plot3(simPowerTest_freqs(i), permute(inactRatioToOrigBySize(i,j,maxInactInd(i,j,countVec),countVec), [4 3 2 1]), log10(simPowerTest_powers(j)), ...
                'x', 'markersize', markerSize, 'color', power_powerCols(j,:))
        end
    end
    xlim([5 12]);  ylim([-10 110])
    plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
%% Calculate predicted inactivation from phenomlogical function

% Set where to make predictions
predictionFreqs = simFreqTest_freqs(1):predFreqStep:simFreqTest_freqs(end);

if ~all(ismember(simFreqTest_freqs, predictionFreqs))
    error('Test freqs missing from prediction calc')
end

if ~all(ismember(simPowerTest_freqs, predictionFreqs))
    error('Test power freqs missing from prediction calc')
end

if ~all(ismember(predictionFreqs, referenceFreqRange))
    error('Test power freqs missing from prediction calc')
end

predictionPowers = simPowerTest_powers;

% Take prediction values from previous
gaussCenterToPred = curveCenter; freqCurve.b1;
gaussMaxToPred = curveMax; freqCurve.a1;
gaussSpreadToPred = curveSpread; freqCurve.c1;

%%% Note, need to get parameters from power from fitted
warning('Power fit not used in pred')

%%% Can allow this to vary if observed does... 
powerThresholdToPred = powerThresholdInterp_ref(ismember(referenceFreqRange, predictionFreqs));

predictedInactivation = zeros(length(predictionFreqs),length(predictionPowers));

%%% Testing using known curves, need to change to fitted.
for i = 1:length(predictionFreqs)
    curveVal = gaussMaxToPred*exp(-(predictionFreqs(i)-gaussCenterToPred).^2/(2*gaussSpreadToPred^2));
    
    for j = 1:length(predictionPowers)
        if predictionPowers(j) > powerThresholdToPred(i)
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

        totalVal = powerVal * curveVal;
        totalVal(totalVal > 100) = 100;     
        totalVal(totalVal < 0) = 0;

        predictedInactivation(i,j) = totalVal;
    end
end

% Remove low lines from power
toRemove = zeros(length(predictionPowers),1,'logical');

for i = 1:length(predictionPowers)
    if all( predictedInactivation(:,i) < predictInactThresh)
        toRemove(i) = 1;
    end
end

powersRemoved = sum(toRemove);

predictionPowers(toRemove) = [];
predictedInactivation(:,toRemove) = [];
referenceInactRatio(:,toRemove,:) = [];

% Remove low lines from freq
toRemove = zeros(length(predictionFreqs),1, 'logical');

for i = 1:length(predictionFreqs)
    if all( predictedInactivation(i,:) < predictInactThresh)
        toRemove(i) = 1;
    end
end

freqsRemoved = sum(toRemove);

predictionFreqs(toRemove) = [];
predictedInactivation(toRemove,:) = [];
referenceInactRatio(toRemove,:,:) = [];

%% Set up for 3D interp

freqsTested = unique([simFreqTest_freqs, simPowerTest_freqs]);

predictionSizes = influenzaSize(sizes2Use);
predictionSizes_dist = influenzaSize_dist(sizes2Use)/sum(influenzaSize_dist)*100;

referenceVol = permute(referenceInactRatio(:,:,sizes2Use), [1 3 2]);
%%% Change to measured data
inactivationVol = referenceVol;

% Clip based on sizes in size map from freq sweep
toRemove = zeros(length(predictionSizes),1,'logical');

for i = 1:length(predictionSizes)
    if all( inactivationVol(:,i,end) == 0) 
        toRemove(i) = 1;
    end
end

sizes2Use(toRemove) = [];
sizes2Remove(toRemove) = [];

predictionSizes(toRemove) = [];
predictionSizes_dist(toRemove) = [];
referenceVol(:,toRemove,:) = [];
inactivationVol(:,toRemove,:) = [];

predictionCountNum = sum(influenzaSize_dist(sizes2Use));
ignoredCountNum = sum(influenzaSize_dist(sizes2Remove));

% Then do on freqs in size map form freq sweep
toRemove = zeros(length(predictionFreqs),1,'logical');

for i = 1:length(predictionFreqs)
    if all( inactivationVol(i,:,end) == 0)
        toRemove(i) = 1;
    end
end

predictionFreqs(toRemove) = [];
predictedInactivation(toRemove,:) = [];
referenceVol(toRemove,:,:) = [];
inactivationVol(toRemove,:,:) = [];

% Get thresholded error for reference
keptInact = size(predictedInactivation);

for i = 1:length(predictionFreqs)
    for j = 1:length(predictionPowers)

        keptInact(i,j) = sum(predictionSizes_dist(:).*referenceVol(i,:,j)'/100);
    end
end

%% Test inactivation

% Test build map using full base and centre line...
[sizeGrid, freqGrid, powerGrid] = meshgrid(1:length(sizes2Use), 1:length(predictionFreqs), 1:length(predictionPowers));

% build exisiting arrays into volumes
measVol = zeros(size(referenceVol));
distVol = zeros(size(referenceVol));

for i = fliplr(1:length(predictionPowers))
    
    tempMeasArray = zeros(size(referenceVol, [1 2]));
    tempInterpArray = zeros(size(referenceVol, [1 2]));
    
    % which freqs are measured
    if i == length(predictionPowers)
        compFreqs = freqsTested;
    else
        compFreqs = simPowerTest_freqs;
    end
    
    freqLinesMeas = [];
    
    % fill out measured array
    for j = 1:length(compFreqs)
        freqInd = find(compFreqs(j) == predictionFreqs);
        
        if ~isempty(freqInd)
            tempMeasArray(freqInd, referenceVol(freqInd,:,i) > interpInactThresh) = -1;
            tempInterpArray(freqInd, :) = -1;
            
            freqLinesMeas = [freqLinesMeas freqInd];
        end
    end    

    measVol(:,:,i) = tempMeasArray;
    
    seedArray =  measVol(:,:,i);
    % If first, do interpolation across
    if i == length(predictionPowers)
        if ~all(ismember(predictionFreqs, freqsTested))
            % first go along each column and fill zeros from top to bottom - makes sure end lines get connected
            for j = 1:length(sizes2Use)
               lineInds = find(seedArray(:,j));

               if ~isempty(lineInds)
                   lineInds = lineInds(1):lineInds(end);

                   seedArray(lineInds(seedArray(lineInds,j) == 0) ,j) = -1;
               end
            end

            % then do interpolant - smooths edges a bit
            [sRef, fRef] = meshgrid(1:length(sizes2Use), 1:length(predictionFreqs));

            indsCheck = find(seedArray == 0 & tempInterpArray == 0);
            indsRef = find(seedArray == -1 | tempInterpArray == -1); %  

            tempInterp = scatteredInterpolant(fRef(indsRef), sRef(indsRef), seedArray(indsRef), 'nearest', 'nearest');

            seedArray(indsCheck) = tempInterp(fRef(indsCheck), sRef(indsCheck));

            seedArray(isnan(seedArray)) = 0;
        end
        
        % this gives full based array
        baseSeedArray = seedArray;
    end
    
    % Do dist map as for 2D
    tempArrayFromLeft = zeros(size(referenceVol, [1 2]));
    tempArrayFromRight = zeros(size(referenceVol, [1 2]));

    levelInds = find(seedArray == -1);
    [levelX, levelY] = ind2sub(size(referenceVol, [1 2]), levelInds);
    
    maxOpenLeft = [];
    maxOpenRight = [];
    
    %%% Right side is not really used anymore
        %%% Include flip before later interp - if right only picked, use that
    
    for j = 1:length(predictionFreqs)
        inds = find(levelX == j);

        if ~isempty(inds)

            yVals = sort(levelY(inds));

            % Check line does not extend across whole map
            if ~(yVals(1) == 1 & yVals(end) == size(referenceVol, 2))
                % Check if first value is not on border
                if yVals(1) > 1
                    % just increases left to right
                    tempArrayFromLeft(levelInds(inds)) = 1:length(inds);
                    
                    maxOpenLeft = [maxOpenLeft length(inds)];
                end

                % Reverse from other side
                %%% Not really needed now
                if yVals(end) < size(referenceVol, 2)
                     %decreases from left to right   
                     tempArrayFromRight(levelInds(inds)) = fliplr(1:length(inds));
                     
                    maxOpenRight = [maxOpenRight length(inds)];
                end

            else
               error('Line extends across whole map') 
            end
        end
    end

    maxOpenLeft = max(maxOpenLeft);
    
    missingLeft = find(tempArrayFromLeft(levelInds) == 0);
    
    % Filling in opposite side
    % Changed to using step back instead of look up table
    flipLeftArray = tempArrayFromLeft;
    for j = 1:length(predictionFreqs)
        inds = find(levelX(missingLeft) == j);
        
        if ~isempty(inds)
            yVals = sort(levelY(missingLeft(inds)));

            % Check line does not extend across whole map
            if yVals(1) ~= 1 | yVals(end) ~= size(referenceVol, 2)
                
                % Step inds down from max
                flipLeftArray(levelInds(missingLeft(inds))) =  (maxOpenLeft-length(inds)+1):maxOpenLeft;
                
                % Removed as in following section, but rarely used?
                % if border is low, step it up
                if flipLeftArray(levelInds(missingLeft(inds(1)))) < 0
                    warning('step up from neg used 1st')
                    flipLeftArray(levelInds(missingLeft(inds))) = flipLeftArray(levelInds(missingLeft(inds))) + 1 - flipLeftArray(levelInds(missingLeft(inds(1))));
                end
            else
               error('Line extends across whole map') 
            end
        end
    end
    
    % If not base, interpolate up as well
    if length(predictionPowers) ~= i & any(flipLeftArray(:) > 0) 
        baseInds = find(baseSeedArray == -1);
        [baseX, baseY] = ind2sub(size(referenceVol, [1 2]), baseInds);
        
        baseArrayFromLeft = distVol(:,:,end);

        % Overlap taken from flipped array
        overlapInds = find(flipLeftArray(baseInds) > 0 & baseArrayFromLeft(baseInds) > 0);
        % Missing filled onto left array so following maxOpenLeft can be updated
        missingLeft = find(tempArrayFromLeft(baseInds) == 0);
        
        % Get low and high overlaps from each row
        referenceLow = zeros(length(predictionFreqs),1);
        referenceHigh = zeros(length(predictionFreqs),1);
        borderLow = zeros(length(predictionFreqs),1);
        borderHigh = zeros(length(predictionFreqs),1);
        
        for j = freqLinesMeas
            inds = find(baseX(overlapInds) == j);
            
            if ~isempty(inds)
                referenceLow(j) = min(baseArrayFromLeft(baseInds(overlapInds(inds))));
                % Flag point on left
                if any(baseY(overlapInds(inds)) == 1)
                    borderLow(j) = 1;
                end
                
                referenceHigh(j) = max(baseArrayFromLeft(baseInds(overlapInds(inds)))); 
                % or high from right borde
                if any(baseY(overlapInds(inds)) == size(referenceVol, 2))
                    borderHigh(j) = 1;
                end
            else
                referenceLow(j) = NaN;
                referenceHigh(j) = NaN;
            end
        end
        
        % find unique rows of base array and give equal freq value in interp 
        [~, inds2Unique, inds2Full] = unique(baseArrayFromLeft,'rows');
        
        tempPredictionFreqs = predictionFreqs;
        for j = 1:length(inds2Full)
            inds = find(inds2Full == inds2Full(j));
            
            tempPredictionFreqs(j) = mean(predictionFreqs(inds));
        end
        
        % Identify if points on border can be used
        notBorderInds = find(referenceLow(freqLinesMeas) > 0 & borderLow(freqLinesMeas) == 0);
        borderInds = find(borderLow(freqLinesMeas));
        if ~isempty(borderInds) & ~isempty(notBorderInds)
            for j = 1:length(borderInds)
                % Using end from not border as left border at high freqs
                    % would be more robust is closest picked
               if referenceLow(freqLinesMeas(borderInds(j))) > referenceLow(freqLinesMeas(notBorderInds(end)))
                   referenceLow(freqLinesMeas(borderInds(j))) = 0;
               end
            end
        end
        
        notBorderInds = find(referenceHigh(freqLinesMeas) > 0 & borderHigh(freqLinesMeas) == 0);
        borderInds = find(borderHigh(freqLinesMeas));
        if ~isempty(borderInds) & ~isempty(notBorderInds)
            for j = 1:length(borderInds)
                % Using 1 from not border as right border at low freqs
               if referenceHigh(freqLinesMeas(borderInds(j))) < referenceHigh(freqLinesMeas(notBorderInds(1)))
                   referenceHigh(freqLinesMeas(borderInds(j))) = 0;
               end
            end
        end
        
        useLow = find(referenceLow(freqLinesMeas) > 0);
        useHigh = find(referenceHigh(freqLinesMeas) > 0);
        
        % deal with missing lines - assumes they are on bounds of good area
        if any(isnan(referenceLow(freqLinesMeas)) | isnan(referenceHigh(freqLinesMeas)))
            blockInds = find(isnan(referenceLow(freqLinesMeas)) | isnan(referenceHigh(freqLinesMeas)));
            useInds = unique([useLow' useHigh']);
            
            freqsToRemove = zeros(length(predictionFreqs),1, 'logical');
            
            for j = 1:length(blockInds)
               if any(useInds > blockInds(j)) & any(useInds < blockInds(j))
                   error('Need to treat')
               elseif any(useInds > blockInds(j))
                   tempTopInds = find(useInds > blockInds(j));
                   tempTopInds = tempTopInds(1);
                   
                   freqsToRemove(1:round((freqLinesMeas(useInds(tempTopInds)) + freqLinesMeas(blockInds(j)))/2)) = 1;
               elseif any(useInds < blockInds(j))
                   tempBotInds = find(useInds < blockInds(j));
                   tempBotInds = tempBotInds(end);
                   
                   freqsToRemove(round((freqLinesMeas(useInds(tempBotInds)) + freqLinesMeas(blockInds(j)))/2):end) = 1;
               end
            end
        else
           freqsToRemove = []; 
        end
        
        % Get unique references to freqs
        [~, tempI, tempBack] = unique([tempPredictionFreqs(freqLinesMeas(useLow))'  referenceLow(freqLinesMeas(useLow))], 'rows','stable');
        useLow = useLow(tempI);
                
        % return different meas results that are on the same compressed line to original freq
        tempPredictionFreqsLow = tempPredictionFreqs;
        for j = 1:length(tempI)
            % If single instance, use original freq
            if sum(tempBack == tempBack(tempI(j))) == 1
                tempPredictionFreqsLow(freqLinesMeas(useLow(j))) =  predictionFreqs(freqLinesMeas(useLow(j)));
            end
        end

        [~, tempI, tempBack] = unique([tempPredictionFreqs(freqLinesMeas(useHigh))' referenceHigh(freqLinesMeas(useHigh))], 'rows','stable');
        useHigh = useHigh(tempI);
        
        % return different meas results that are on the same compressed line to original freq
        tempPredictionFreqsHigh = tempPredictionFreqs;
        for j = 1:length(tempI)
            % If single instance, use original freq
            if sum(tempBack == tempBack(tempI(j))) == 1
                tempPredictionFreqsHigh(freqLinesMeas(useHigh(j))) =  predictionFreqs(freqLinesMeas(useHigh(j)));
            end
        end
        
%         [tempPredictionFreqsLow' tempPredictionFreqs' tempPredictionFreqsHigh']
        
        % interp internal with linear and extrap from border with nearest
        if length(useLow) > 1
            referenceLow(freqLinesMeas(useLow(1)):freqLinesMeas(useLow(end))) = round(interp1(tempPredictionFreqsLow(freqLinesMeas(useLow)), referenceLow(freqLinesMeas(useLow)),...
                tempPredictionFreqsLow(freqLinesMeas(useLow(1)):freqLinesMeas(useLow(end))), distInterpType));
            referenceLow(1:freqLinesMeas(useLow(1))) = referenceLow(freqLinesMeas(useLow(1)));
            referenceLow(freqLinesMeas(useLow(end)):end) = referenceLow(freqLinesMeas(useLow(end)));
        else
            referenceLow(:) = referenceLow(freqLinesMeas(useLow));
        end
        
        if length(useHigh) > 1
            referenceHigh(freqLinesMeas(useHigh(1)):freqLinesMeas(useHigh(end))) = round(interp1(tempPredictionFreqsHigh(freqLinesMeas(useHigh)), referenceHigh(freqLinesMeas(useHigh)),...
                tempPredictionFreqsHigh(freqLinesMeas(useHigh(1)):freqLinesMeas(useHigh(end))), distInterpType));
            referenceHigh(1:freqLinesMeas(useHigh(1))) = referenceHigh(freqLinesMeas(useHigh(1)));
            referenceHigh(freqLinesMeas(useHigh(end)):end) = referenceHigh(freqLinesMeas(useHigh(end)));
        else
            referenceHigh(:) = referenceHigh(freqLinesMeas(useHigh));
        end

        % Remove inds for missing lines
        if ~isempty(freqsToRemove)
            referenceLow(freqsToRemove) = 0;
            referenceHigh(freqsToRemove) = 0;
        end
        
        for j = 1:length(predictionFreqs)
            if all(tempArrayFromLeft(j,:) == 0)
                inds = find(baseX(missingLeft) == j & baseArrayFromLeft(baseInds(missingLeft)) >= referenceLow(j) & ...
                    baseArrayFromLeft(baseInds(missingLeft)) <= referenceHigh(j));

                if ~isempty(inds)
                    yVals = sort(baseY(missingLeft(inds)));

                    % Check line does not extend across whole map
                    if yVals(1) ~= 1 | yVals(end) ~= size(referenceVol, 2)

                        % Check if first value is not on left border
                        if yVals(1) > 1
                            % just count up normally
                            tempArrayFromLeft(baseInds(missingLeft(inds))) = 1:length(inds);
                            
                            if maxOpenLeft < length(inds)
                                maxOpenLeft = length(inds);
                            end
                        else
                             % Step inds down from max
                             tempArrayFromLeft(baseInds(missingLeft(inds))) =  (maxOpenLeft-length(inds)+1):maxOpenLeft;

                            % Removed as this seemed to skew the interp map left
                            % if border is low, step it up
                            if tempArrayFromLeft(baseInds(missingLeft(inds(1)))) < 0
                                warning('step up from neg used 2nd')
                                tempArrayFromLeft(baseInds(missingLeft(inds))) = tempArrayFromLeft(baseInds(missingLeft(inds))) + 1 - tempArrayFromLeft(baseInds(missingLeft(inds(1))));
                            end
                        end
                    else
                        error('Line extends across whole map') 
                    end
                end
            end
        end
    else
       tempArrayFromLeft = flipLeftArray; 
    end
    
    distVol(:,:,i) = tempArrayFromLeft; %distanceArrayFromLeft;
end

%%% Switch to using inactivated volume below
measInds = find(measVol == -1);

allInds = find(distVol);

length(find(distVol == 0 & referenceVol > interpInactThresh))

interpVol = zeros(size(referenceVol));

% Now do interp
if length(simPowerTest_freqs) > 1
    % nearest best for both
    inactInterp = scatteredInterpolant(distVol(measInds), freqGrid(measInds), powerGrid(measInds), referenceVol(measInds), 'nearest', 'nearest');

    interpVol(allInds) = inactInterp(distVol(allInds), freqGrid(allInds), powerGrid(allInds));
else
    % triangulation fails x included - nearest best
    inactInterp = scatteredInterpolant(distVol(measInds), powerGrid(measInds), referenceVol(measInds), 'nearest', 'nearest');

    interpVol(allInds) = inactInterp(distVol(allInds), powerGrid(allInds));
end

% above will overwrite measured inds if powers didn't match, so add back
interpVol(measInds) = referenceVol(measInds);

interpVol(interpVol > 100) = 100;
interpVol(interpVol < 0) = 0;
interpVol(isnan(interpVol(:))) = 0;

figure;
for i = 1:length(predictionPowers)
    subplot(length(predictionPowers),4,i*4-3)
    imshow(referenceVol(:,:,i)/100)
    
    subplot(length(predictionPowers),4,i*4-2)
    imshow(interpVol(:,:,i)/100)
    
    subplot(length(predictionPowers),4,i*4-1)
    imshow(distVol(:,:,i)/max(max(distVol(:,:,i))))
    
    subplot(length(predictionPowers),4,i*4)
    imshow(measVol(:,:,i)*-1)
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
imshow(keptInact/100)
title('Ref kept inact')

subplot(4,4,5)
imshow(interpInactArray/100)
title('From volume')

subplot(4,4,6)
imshow(abs(predictedInactivation-interpInactArray)/10)
title('Dif full to vol')

subplot(4,4,7)
imshow(abs(keptInact-interpInactArray)/10)
title('Dif kept to vol')

subplot(4,4,8); hold on
plot(predictionFreqs, sum(abs(predictedInactivation-interpInactArray),2),'g')
title('error vs full across')

subplot(4,4,12); hold on
plot(predictionFreqs, sum(abs(keptInact-interpInactArray),2),'g')
title('error vs kept across')

subplot(4,4,10); hold on
plot(predictionPowers, sum(abs(predictedInactivation-interpInactArray),1),'g')
title('error vs full down')

subplot(4,4,11); hold on
plot(predictionPowers, sum(abs(keptInact-interpInactArray),1),'g')
title('error vs kept down')

sseVol_full = sum((predictedInactivation(:)-interpInactArray(:)).^2);
sseVol_thresh = sum((keptInact(:)-interpInactArray(:)).^2);

[sseVol_full, sseVol_thresh]