compare_influenza_data
close all
%% - generate parameters for size - would be from SEM

influenzaSize_mean = 100*10^-9;
influenzaSize_std = 20*10^-9/3; % Assume limits are at 3 standard deviations

initialCountNum = 5000;

influenzaSize_samples = randn(initialCountNum,1)*influenzaSize_std + influenzaSize_mean;

% How accurately can size be determined from TEM? surely not < 1 nm... 

distStep = 2; % Also percentage as across 100
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
    % Later - need to add replication and calibration test on distribution

absStd= 5; % std on absolute, not relative values (lower SNR on low inactiviation)

% phenomenlogical parameters
    % For freq
        curveMax = 100;
        curveCenter = 8.5; 
        curveSpread = 2.3;

    % For power - single threhsold
        % powerThreshold = 45; %45
        % For power - linear interpolation between points then constant
        powerThreshold = [35 55]; 
        powerThresholdFreqs = [8 9]

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
simFreqTest_freqs = 1:1:20;
simFreqTest_power = 630;
simFreqTest_time = 1000;

% Phase 1.2 - identify power response (minimum effective) at center freq
powerLogRange = [1:0.25:2.25]; % coarse search
%%% Note neccersary to duplicate powers from previous step i.e. start at 1.55 and 1.7
% powerLogRange = 1.5:0.05:1.75; % For fine search. works on 45 or 33
simPowerTest_powers = round(10.^(powerLogRange))
simPowerTest_time = simFreqTest_time;

%%% Should do second step to get fine values

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

% Make simulated plots
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

    totalVal = powerVal * curveVal;
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

simPowerTest_freqs = freqCurve.b1;

% Look at distributions - take for first result
figure; hold on

%%% Make a function for this, will use later
for i = 1:length(simFreqTest_freqs)
    
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

        subplot(2,10,i); hold on
        plot(influenzaSize*10^9, activeDistDist+inactiveDist, 'r')

        plot(influenzaSize*10^9, activeDistDist, 'k')
    end
end

figure;
subplot(1,3,1); hold on
plot(dataInactFreq2016(:,1), dataInactFreq2016(:,2), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
plot(simFreqTest_FreqRef(:), simFreqTest_inact(:), 'rd', 'markersize', 4)
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

plot(freqRef, powerInd, 'rd', 'markersize', 8)

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
        
        totalVal = powerVal * curveVal;
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
subplot(3,length(simPowerTest_powers),1); hold on

for i = 1:length(simPowerTest_freqs)
    for j = 1:length(simPowerTest_powers)
        toPlot = setxor(countVec, 1:nReps);

        plot(log10(permute(simPowerTest_powerRef(i,j,toPlot), [3 2 1])), permute(simPowerTest_inact(i,j,toPlot), [3 2 1]), 'ro', 'markersize', markerSize)

        plot(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(simPowerTest_inact(i,j,countVec), [3 2 1]), 'r*', 'markersize', markerSize)
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

plot(log10(permute(simPowerTest_powerRef(1,:,1), [3 2 1])), permute(simPowerTest_inactRef(1,:,1), [3 2 1]), 'm', 'linewidth', 2)



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
for i = 1
    for j = fliplr(1:length(simPowerTest_powers))

        subplot(3,length(simPowerTest_powers),length(simPowerTest_powers)+j); hold on
        plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-b', 'linewidth', 2)

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
            subplot(3,length(simPowerTest_powers),length(simPowerTest_powers)+j); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'k')

            % inactivation spectrum
            subplot(3,length(simPowerTest_powers),2*length(simPowerTest_powers)+j); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,:), [4 3 2 1]), 'k');
            ylim([-50 100])
        end
    end
end

% inact from counted intact summed
subplot(3,length(simPowerTest_powers),4); hold on
for i = 1
    for j = 1:length(simPowerTest_powers)
        plot(log10(permute(simPowerTest_powerRef(1,j,countVec), [3 2 1])), permute(totalInact(1,j,countVec), [4 3 2 1]), 'rx', 'markersize', markerSize)
    end
end

for i = 1
    for j = 1:length(simPowerTest_powers)
        plot(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(simPowerTest_inact(i,j,countVec), [3 2 1]), 'r*', 'markersize', markerSize)
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])


% Inact from counted both across range
subplot(3,length(simPowerTest_powers),2); hold on

for i = 1
    for j = 1:length(simPowerTest_powers)
        plot(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4, 3 2 1]),...
            'rx', 'markersize', markerSize)
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Inact from counted intact across range
subplot(3,length(simPowerTest_powers),3); hold on
for i = 1
    for j = 1:length(simPowerTest_powers)
        plot(log10(permute(simPowerTest_powerRef(i,j,countVec), [3 2 1])), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]),...
            'rx', 'markersize', markerSize)
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Plot inactivation width

%%% Could use thickness or color to indicate number inactivated
subplot(3,length(simPowerTest_powers),5); hold on
for i = 1
    for j = 1:length(simPowerTest_powers)
        for k = countVec
            if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2)
            line(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9,...
                'color', 'r', 'linewidth', 2)
            else
                if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100
                    plot(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9,...
                        'r.', 'markersize', 8)
                else
                    plot(log10(permute(simPowerTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9,...
                        'ro', 'markersize', 4)
                end
            end
        end
    end
end
ylim([80 120])
xlim([1 3])
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
        
        totalVal = powerVal * curveVal;
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
subplot(3,length(simFineTest_freqs),1); hold on

for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)
        toPlot = setxor(countVec, 1:nReps);
        plot3(log10(permute(simFineTest_powerRef(i,j,toPlot), [3 2 1])), permute(simFineTest_inact(i,j,toPlot), [3 2 1]), permute(simFineTest_freqRef(i,j,toPlot), [3 2 1]), ...
            'o', 'markersize', markerSize, 'color', fine_freqCols(i,:));
        
        plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(simFineTest_inact(i,j,countVec), [3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
            '*', 'markersize', markerSize, 'color', fine_freqCols(i,:));
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% plot(log10(permute(simFineTest_powerRef(1,:,1), [3 2 1])), permute(simFineTest_inactRef(1,:,1), [3 2 1]), 'm', 'linewidth', 2)

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
    
    subplot(3,length(simFineTest_freqs),length(simFineTest_freqs)+i); hold on
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
            subplot(3,length(simFineTest_freqs),length(simFineTest_freqs)+i); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'color', fine_powerCols(j,:))

            % inactivation spectrum
            subplot(3,length(simFineTest_freqs),2*length(simFineTest_freqs)+i); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,:), [4 3 2 1]), 'color', fine_powerCols(j,:));
            ylim([-50 100])
        end
    end
end

% inact from counted intact summed
subplot(3,length(simFineTest_freqs),4); hold on
for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)
        plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(totalInact(i,j,countVec), [4 3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
            'x', 'markersize', markerSize, 'color', fine_freqCols(i,:))
    end
end

for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)       
        plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(simFineTest_inact(i,j,countVec), [3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
            '*', 'markersize', markerSize, 'color', fine_freqCols(i,:));
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Inact from counted both across range
subplot(3,length(simFineTest_freqs),2); hold on

for i = 1:length(simFineTest_freqs)
    for j = 1:length(simFineTest_powers)
        plot3(log10(permute(simFineTest_powerRef(i,j,countVec), [3 2 1])), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4 3 2 1]), permute(simFineTest_freqRef(i,j,countVec), [3 2 1]), ...
            'x', 'markersize', markerSize, 'color', fine_freqCols(i,:))
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Inact from counted intact across range
subplot(3,length(simFineTest_freqs),3); hold on

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
subplot(3,length(simFineTest_freqs),5); hold on
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
                else
                    plot3(log10(permute(simFineTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9, permute(simFineTest_freqRef(i,j,countVec), [3 2 1]),...
                        'o', 'color', fine_freqCols(i,:), 'markersize', 4)
                end
            end
        end
    end
end

ylim([80 120])
xlim([1 3])

%% Phase 1.4 - test time

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

for i = 1:length(simTimeTest_freqs)

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
subplot(4,length(simTimeTest_powers),1); hold on

for i = 1:length(simTimeTest_freqs)
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
subplot(4,length(simTimeTest_powers),length(simTimeTest_powers)+1); hold on

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

inactRatioBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount, length(influenzaSize));

inactFractionBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount, length(influenzaSize));

activeBySize = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount, length(influenzaSize));

totalInact = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount);


inactInds = cell(length(simTimeTest_freqs), length(simTimeTest_powers));

maxInactInd = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount); % if more than 1 equal, will just be one 

maxInactBounds = zeros(length(simTimeTest_freqs), length(simTimeTest_powers), nrepsCount, 2);

%%% Make a function for this, will use later
for i = 1
    for j = fliplr(1:length(simTimeTest_powers))

        subplot(3,length(simTimeTest_powers),length(simTimeTest_powers)+j); hold on
        plot(influenzaSize*10^9, influenzaSize_dist/initialCountNum*100, '-b', 'linewidth', 2)

        for k = countVec
            tempSize_samples = randn(testCountNum,1)*influenzaSize_std + influenzaSize_mean;

            tempSize_samples = sort(tempSize_samples, 'descend');

            tempSize_freqs = 1./(tempSize_samples/2)*resSlope;

            samplesIntact = ones(testCountNum, 1);

            % Start with reference before we figure out noise... 
            numToInact = round(testCountNum*simTimeTest_inactRef(i,j,k)/100);

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

            inactRatioBySize(i,j,k,:) = (1 - (activeDistDist/testCountNum)./(influenzaSize_dist/initialCountNum))*100;

            inactRatioBySize(i,j,k,zeroRatioInds) = 0;

            inactFractionBySize(i,j,k,:) = (1 - (activeDistDist./(activeDistDist+inactiveDist)))*100;

            inactFractionBySize(i,j,k,zeroFractionInds) = 0;

            activeBySize(i,j,k,:) = activeDistDist/testCountNum*100;

            totalInact(i,j,k) = (1-sum(activeDistDist)/testCountNum)*100;

            % get inds and ranges on inactivation spectrum
            
            % find max
            if length(simTimeTest_powers) == j
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
            subplot(3,length(simTimeTest_powers),length(simTimeTest_powers)+j); hold on

            plot(influenzaSize*10^9, activeDistDist/testCountNum*100, 'k')

            % inactivation spectrum
            subplot(3,length(simTimeTest_powers),2*length(simTimeTest_powers)+j); hold on

            plot(influenzaSize*10^9, permute(inactRatioBySize(i,j,k,:), [4 3 2 1]), 'k');
            ylim([-50 100])
        end
    end
end

% inact from counted intact summed
subplot(3,length(simTimeTest_powers),4); hold on
for i = 1
    for j = 1:length(simTimeTest_powers)
        plot(log10(permute(simTimeTest_powerRef(1,j,countVec), [3 2 1])), permute(totalInact(1,j,countVec), [4 3 2 1]), 'rx', 'markersize', markerSize)
    end
end

for i = 1
    for j = 1:length(simTimeTest_powers)
        plot(log10(permute(simTimeTest_powerRef(i,j,countVec), [3 2 1])), permute(simTimeTest_inact(i,j,countVec), [3 2 1]), 'r*', 'markersize', markerSize)
    end
end

plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])


% Inact from counted both across range
subplot(3,length(simTimeTest_powers),2); hold on

for i = 1
    for j = 1:length(simTimeTest_powers)
        plot(log10(permute(simTimeTest_powerRef(i,j,countVec), [3 2 1])), permute(inactFractionBySize(i,j,countVec, maxInactInd(i,j,countVec)), [4, 3 2 1]),...
            'rx', 'markersize', markerSize)
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Inact from counted intact across range
subplot(3,length(simTimeTest_powers),3); hold on
for i = 1
    for j = 1:length(simTimeTest_powers)
        plot(log10(permute(simTimeTest_powerRef(i,j,countVec), [3 2 1])), permute(inactRatioBySize(i,j,countVec,maxInactInd(i,j,countVec)), [4 3 2 1]),...
            'rx', 'markersize', markerSize)
    end
end
plot(log10(dataInactPower2016(:,1)), dataInactPower2016(:,2), '-', 'linewidth',2, 'color', [0.5 0.5 0.5])
xlim([1 3])

% Plot inactivation width

%%% Could use thickness or color to indicate number inactivated
subplot(3,length(simTimeTest_powers),5); hold on
for i = 1
    for j = 1:length(simTimeTest_powers)
        for k = countVec
            if maxInactBounds(i,j,k,1) ~= maxInactBounds(i,j,k,2)
            line(log10(permute(simTimeTest_powerRef(i,j,k), [3 2 1]))*[1 1], permute(maxInactBounds(i,j,k,:), [4 3 2 1])*10^9,...
                'color', 'r', 'linewidth', 2)
            else
                if inactRatioBySize(i,j,k,maxInactInd(i,j,k)) == 100
                    plot(log10(permute(simTimeTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9,...
                        'r.', 'markersize', 8)
                else
                    plot(log10(permute(simTimeTest_powerRef(i,j,k), [3 2 1])), permute(maxInactBounds(i,j,k,1), [4 3 2 1])*10^9,...
                        'ro', 'markersize', 4)
                end
            end
        end
    end
end
ylim([80 120])
xlim([1 3])