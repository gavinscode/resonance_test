clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [4 6]; [1 4 6];    

nSize = length(sizesToUse);

simWithInterpSizeDist = 1; % If not set, uses original distribution for simulation

varyCharge = 0; % in simulation, using quadratic model?

sizeSteps = 1/10^9; % in m

fitRealData = 1; % if not set, fits simulated data

numberOfModes = 3;

% fitWithSizeType = 'Single';
fitWithSizeType = 'Original';
% fitWithSizeType = 'FullInterp';

typeOfFit = '3Coeffs';

% typeOfFit = 'IntensityPerSize';
    constrainToIncreasing = 0; %fmincon will be used with linear inequality

typeOfFit = 'IntensityFunction';
%     typeOfFunction = 'Linear';
    typeOfFunction = 'Quadratic';
%     typeOfFunction = 'Cubic';

useHandFreqLimits = 0;

if ~useHandFreqLimits
    frequencyRange_rad = (100:5:400)*10^9*2*pi;
end

frequencyRangePlot_rad = (100:1:400)*10^9*2*pi;

% First do for each indvidually and store data as we go
diameterDistGroup = cell(nSize,1);
numberGroup = cell(nSize,1);

bulkQFGroup = zeros(nSize,1);
bulkIntensityGroup = zeros(nSize,1);
bulkResonantSlopeGroup = zeros(nSize,1);

% will only do group fitting on three param options.
parametersToUseGroup = zeros(nSize,3);

dataToFitGroup = cell(nSize,1);
frequencyGroup = cell(nSize,1);

for iSize = 1:nSize
    sizeIndex = sizesToUse(iSize);

    % Get size dist
    diameterDistOrig = nanocrystalSizeDist{sizeIndex, 1};
    
    countDistOrig = nanocrystalSizeDist{sizeIndex, 2};
    
    sizeFrequencyOrig = countDistOrig/sum(countDistOrig);
    
    [countDistInterp, diameterDistInterp] = interpolatescaleddistribution(countDistOrig, diameterDistOrig, sizeSteps);
    
    sizeFrequencyInterp = countDistInterp/sum(countDistInterp);
    
    if simWithInterpSizeDist
        diameterDistSim = diameterDistInterp;

        countDistSim = countDistInterp;
        
    else
        diameterDistSim = diameterDistOrig;

        countDistSim = countDistOrig;
    end

    sizeFrequencySim = countDistSim/sum(countDistSim);
    
    if useHandFreqLimits
        frequencyRange_rad = (lowerFreqToFit_1st(sizeIndex):5:upperFreqToFit_1st(sizeIndex))*10^9*2*pi;
    end
    
    measuredExtinctionCurve = ExCrossSecCurve_m2{sizeIndex};
    % Interpolate to currenct frequncies
    measuredExtinctionCurve = interp1(curveFrequncy*2*pi, measuredExtinctionCurve, frequencyRange_rad, 'linear', 'extrap');
    
    measuredExtinctionCurve(measuredExtinctionCurve < 0) = 0;
    
    % Calculate parameters from bulk
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;
    
    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;
    
    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
    
    reducedMass = coreMass*shellMass/(coreMass + shellMass);
    
    %%% For now, use Q and q from bulk
    resonance_rad = nanocrystalFreqResonance_hz(sizeIndex)*2*pi;
    
    % Q is an unknown, unsure how it should vary with size
    bulkQF = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);

    systemSpring = resonance_rad^2*reducedMass;

    systemDamp = resonance_rad*reducedMass/bulkQF;

    % q is an unknown and should vary with the square of size 
    bulkCharge = sqrt(nanocrystalThetaEx_m2(sizeIndex) * resonance_rad * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED/bulkQF);  

    % Get slope of main size and 1st mode
    resonantSlope = resonance_rad/(1/(nanocrystalSize_m(sizeIndex)/2));
    
    % Take multipliers of other modes from 1st
    slopeMult = [nanocrystal2ndFreqResonance_hz(sizeIndex)/nanocrystalFreqResonance_hz(sizeIndex), ...
        nanocrystal3rdFreqResonance_hz(sizeIndex)/nanocrystalFreqResonance_hz(sizeIndex)];
    
    slopeMult(isnan(slopeMult)) = [];
    
    % Absorbtion and extinction for bulk    
    [bulkAbsorbtion, bulkExtinctionCrossSection, bulkPower] = calculatesphereabsorbtion(frequencyRangePlot_rad, resonance_rad, ...
        reducedMass, bulkQF, bulkCharge, nanocrystalNumber(sizeIndex), apertureArea, 1, []);
    
    [~, testExtinction] = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, resonance_rad, ...
        bulkCharge/sqrt(reducedMass), bulkQF, nanocrystalNumber(sizeIndex), []);
    
    % Store for group fitting
    bulkQFGroup(iSize) = bulkQF;
    bulkIntensityGroup(iSize) = bulkCharge/sqrt(reducedMass);
    bulkResonantSlopeGroup(iSize) = resonantSlope;
    
    % Make coarse synthetic and plot against original and single
    %%% Fixed Q and q - should be easy to fit...
    
    resonanceBySize_rad = zeros(length(diameterDistSim), 1);
    
    reducedMassBySize = zeros(length(diameterDistSim), 1);
    
    if varyCharge
        % Set charge in e as quadractic fn of diameter in nm
        chargeToUse = 2*(diameterDistSim*10^9).^2*(1.602176634*10^-19);
    else
        chargeToUse = bulkCharge;
    end
    
    for jDiameter = 1:length(diameterDistSim)

        coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDistSim(jDiameter)/ ...
            nanocrystalSize_m(sizeIndex);
        
        coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

        totalVolume = 4/3*pi*(diameterDistSim(jDiameter)/2).^3;

        coreMass = coreVolume*CdSeDensity_kgpm3;

        shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
        
        % This will go to zero is shell has zero mass - no absorbtion occurs
        reducedMassBySize(jDiameter) = coreMass*shellMass/(coreMass + shellMass);
        
        % Get resonant frequncy
        resonanceBySize_rad(jDiameter) = 1./(diameterDistSim(jDiameter)/2) * resonantSlope;
    end
    
%     figure; subplot(1,3,1);
%     plot(diameterDistSim*10^9, chargeToUse./sqrt(reducedMassBySize)')
%     subplot(1,3,2);
%     plot(diameterDistSim*10^9, chargeToUse'./sqrt(reducedMassBySize).*sizeFrequency')
%     subplot(1,3,3);
%     plot(diameterDistSim*10^9, reducedMassBySize);
    
    % Calculate extinction for sphere mixture
    [~, simulatedExtinctionCrossSection] = calculatespheremixtureabsorbtion(frequencyRange_rad, resonanceBySize_rad, ...
        reducedMassBySize, bulkQF, chargeToUse, nanocrystalNumber(sizeIndex)*sizeFrequencySim, apertureArea, 1, []);
    
    figure;
    plot(frequencyRange_rad/2/pi/10^9, measuredExtinctionCurve/10^-21, 'm-'); hold on
    
    plot(frequencyRangePlot_rad/2/pi/10^9, bulkExtinctionCrossSection/10^-21, 'g')

    plot(frequencyRangePlot_rad/2/pi/10^9, testExtinction/10^-21, 'r--')
    
    plot(frequencyRange_rad/2/pi/10^9, simulatedExtinctionCrossSection/10^-21, 'g')
    
    % Choose data to fit
    if fitRealData
        dataToFit = measuredExtinctionCurve;
    else
        dataToFit = simulatedExtinctionCrossSection;
    end
    
    % Dimensions of dataToFit must match calculateresonantormixtureabsorbtion output
    % otherwise solver fails quietly
    if size(dataToFit,2) > size(dataToFit, 1)
       dataToFit = dataToFit'; 
    end
    
    % Choose sizes and frequency to use in fit
    switch fitWithSizeType
        case 'Single'
            diametersToUse = nanocrystalSize_m(sizeIndex);
            numberToUse = nanocrystalNumber(sizeIndex);
            
        case 'Original'
            diametersToUse = diameterDistOrig;
            numberToUse = nanocrystalNumber(sizeIndex)*sizeFrequencyOrig;
            
        case 'FullInterp'
            diametersToUse = diameterDistInterp;
            numberToUse = nanocrystalNumber(sizeIndex)*sizeFrequencyInterp;
    end
    
    % Store for group fitting
    diameterDistGroup{iSize} = diametersToUse;
    numberGroup{iSize} = numberToUse;
    
    dataToFitGroup{iSize} = dataToFit;
    frequencyGroup{iSize} = frequencyRange_rad;
    
    % Limit to number of modes with indicated peaks, otherwise top mode floats
    if length(slopeMult) + 1 < numberOfModes
        modesToTest = length(slopeMult) + 1;
    else
        modesToTest = numberOfModes;
    end
    
    switch modesToTest
        
        case 1
            slopeMultX0 = [];
            slopeMultLb = [];
            slopeMultUb = [];
            
            slopeRef = [];
            
        case 2
            %%% Using limits from slop map tester...
            slopeMultX0 = 1.6;
            slopeMultLb = 1.4;
            slopeMultUb = 2.0;
            
            slopeRef = [2];
        case 3
            slopeMultX0 = [1.6 2.8];
            slopeMultLb = [1.4 2.3];
            slopeMultUb = [2.0 3.1];
            
            slopeRef = [2 2];
    end
    
    % Set up fit
    switch typeOfFit
    
        % Use solver for 3 coefficients
        case '3Coeffs'
            % Parameters rescalled to similar ranges
            x0 = [resonantSlope/1000, slopeMultX0, bulkCharge/sqrt(reducedMass)*10^6, bulkQF*ones(1, modesToTest)];
                  
            lb = [resonantSlope/1000/2, slopeMultLb, 0, 0*ones(1, modesToTest)];
            ub = [resonantSlope/1000*2, slopeMultUb, 20, 20*ones(1, modesToTest)];
            
            positionRef = [1 slopeRef 3 4*ones(1, modesToTest)];
            
            f = @(x)(calculateerror_equalcoefficients(x, dataToFit, diametersToUse, ...
                numberToUse, frequencyRange_rad, positionRef));
    
        % Use solver allowing for varying intensity    
        case 'IntensityPerSize'
            warning('Multi mode not implemented')
            
            %%% Soln values quite sensitive to initial parameter choices
                % Acceptable space deffinetely worth considering...
            x0 = [resonantSlope/1000, bulkCharge/sqrt(reducedMass)*10^6*ones(1,length(diametersToUse)),...
                bulkQF];

            lb = [resonantSlope/1000/2 zeros(1,length(diametersToUse)) 1];
            ub = [resonantSlope/1000*2 20*ones(1,length(diametersToUse))  20];

            % Set up linear constraints on indvidual charge
            A = zeros(length(x0), length(x0));

            % Each should be large than precedding
            for jCharge = 3:length(x0)-1
               A(jCharge, jCharge-1:jCharge) = [1 -1]; 
            end

            %%% Would also be reasonable to let all larger than a given size go
            %%% to zero (indicating none of that size.)

            b = zeros(length(x0), 1);

            f = @(x)(calculateerror_varyingintensity(x, dataToFit, diametersToUse, ...
                numberToUse, frequencyRange_rad, constrainToIncreasing));
            
        case 'IntensityFunction'
%             x0 = [resonantSlope/1000, 0.01, bulkQF];
% 
%             lb = [resonantSlope/1000/2 0.001 1];
%             ub = [resonantSlope/1000*2 0.5 20];
            
            x0 = [resonantSlope/1000, slopeMultX0, 0.01, bulkQF*ones(1, modesToTest)];
                  
            lb = [resonantSlope/1000/2, slopeMultLb, 0.001, ones(1, modesToTest)];
            ub = [resonantSlope/1000*2, slopeMultUb, 0.5, 20*ones(1, modesToTest)];
            
            positionRef = [1 slopeRef 3 4*ones(1, modesToTest)];
            
            f = @(x)(calculateerror_intensityfunction(x, dataToFit, diametersToUse, ...
                numberToUse, frequencyRange_rad, typeOfFunction, positionRef));
    end

    
    % Use lower bound at zero for all
    if constrainToIncreasing
        solution = fmincon(f, x0, A, b, [], [], lb, ub)
        
        (A*solution')'
    else
        options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10, ...
            'MaxFunctionEvaluations', 10^6, 'MaxIterations', 10^4);
        
        solution = lsqnonlin(f, x0, lb, ub, options)
    end
    
    solution - x0
    
    %%% Need to remake...
%     if length(solution) == 3
%         parametersToUseGroup(iSize, :) = solution;
%     end
    
    % Unpack solution
    slopeInds = find(positionRef == 1);
    slopeMultInds = find(positionRef == 2);
    chargeInds = find(positionRef == 3);
    QFInds = find(positionRef == 4); 
    
    solvedSlope = solution(slopeInds)*1000;
    
    solvedSlopeMult = solution(slopeMultInds);
    
    solvedCharges = solution(chargeInds)/10^6;
    
    solvedQFs = solution(QFInds);
    
    % Calculate parameters solved for
    resonance_rad = 1./(diametersToUse/2) * solvedSlope;
    resonanceBase_rad = resonance_rad;
    
    % Intensity depends on type of solution
    switch typeOfFit
        case '3Coeffs'
            intensityValues = solvedCharges;
            
        case 'IntensityPerSize'
            intensityValues = solvedCharges;
            
        case 'IntensityFunction'
            
            switch typeOfFunction
                case 'Linear'
                    % Diameter in nm
                    refSize = diametersToUse*10^9;

                case 'Quadratic'
                    % Area in nm^2
                    refSize =  4*pi*(diametersToUse/2).^2*10^18;

                case 'Cubic'
                    %volume in nm^3
                    refSize = 4/3*pi*(diametersToUse/2).^3*10^27;
            end
            
            intensityValues = refSize*solvedCharges;
    end
    
    % Check predicted cross section for each mode
    
    abosrbtion = zeros(length(frequencyRangePlot_rad),1);
    
    % Do for each mode and combine absorbtions
    for iMode = 1:modesToTest
        if iMode > 1
            resonance_rad = resonanceBase_rad*solvedSlopeMult(iMode-1);
        end
        
        tempAbs = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, resonance_rad, ...
            intensityValues, solvedQFs(iMode), numberToUse, []);
        
        abosrbtion = 1 - (1 - abosrbtion) .* (1 - tempAbs);
    end

    extinctionCrossSection = -1/sum(numberToUse).*log(1-abosrbtion);
    
    plot(frequencyRangePlot_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'r--')
end

%% Run on combined!
x0 = mean(parametersToUseGroup);

lb = [x0(1)/2 0.001 1];
ub = [x0(1)*2 0.5 20];

f = @(x)(calculateerror_intensityfunction(x, dataToFitGroup, diameterDistGroup, ...
    numberGroup, frequencyGroup, 10^21, typeOfFunction));

options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 10^6, 'MaxIterations', 10^4);

solution = lsqnonlin(f, x0, lb, ub, options)

solution - x0
    
parametersToUse = solution;

solution = [solution(1)*1000 solution(2:end-1)/10^6 solution(end)];

figure;

mode1Data = cell(nSize,1);

%Plot for each size
for iSize = 1:nSize
    subplot(2,nSize,iSize); hold on
    
    diametersToUse = diameterDistGroup{iSize};

    resonance_rad = 1./(diametersToUse/2) * solution(1);
    
    switch typeOfFit
        case '3Coeffs'
            intensityValues = solution(2);
            
        case 'IntensityPerSize'
            intensityValues = solution(2:end-1);
            
        case 'IntensityFunction'
            
            switch typeOfFunction
                case 'Linear'
                    % Diameter in nm
                    refSize = diametersToUse*10^9;

                case 'Quadratic'
                    % Area in nm^2
                    refSize =  4*pi*(diametersToUse/2).^2*10^18;

                case 'Cubic'
                    %volume in nm^3
                    refSize = 4/3*pi*(diametersToUse/2).^3*10^27;
            end
            
            intensityValues = refSize*solution(2);
    end
    
    % Check predicted cross section
    [absorbtion, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, resonance_rad, ...
        intensityValues, solution(end), numberGroup{iSize}, []);
    
    mode1Data{iSize} = extinctionCrossSection;
    
    plot(frequencyRangePlot_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'r--')
    
    sizeIndex = sizesToUse(iSize);
    
    % Get full extinction curve
    measuredExtinctionCurve = interp1(curveFrequncy*2*pi, ExCrossSecCurve_m2{sizeIndex}, frequencyRangePlot_rad, 'linear', 'extrap');
    
    plot(frequencyRangePlot_rad/2/pi/10^9, measuredExtinctionCurve/10^-21, 'm-'); hold on
    
    % Forgot to take difference of extinction properly...
    measAbsCurve = (1-exp(-measuredExtinctionCurve*nanocrystalNumber(sizeIndex)/1));
    
    diffAbsCurve = 1 - (1 - measAbsCurve)./(1 - absorbtion');
    
    diffExtCurve = -1/nanocrystalNumber(sizeIndex).*log(1-diffAbsCurve);
    
    % Difference curve
    plot(frequencyRangePlot_rad/2/pi/10^9, (diffExtCurve)/10^-21, 'b-'); hold on
    
    % Set up offset curve and frequncies of second peak
    if useHandFreqLimits
        frequencyRange_rad = (lowerFreqToFit_2nd(sizeIndex):5:upperFreqToFit_2nd(sizeIndex))*10^9*2*pi;
    end
    
    [tempAbsCalc] = calculateresonantormixtureabsorbtion(frequencyRange_rad, resonance_rad, ...
        intensityValues, solution(end), numberGroup{iSize}, []);
    
    tempMeasEx = interp1(curveFrequncy*2*pi, ExCrossSecCurve_m2{sizeIndex}, frequencyRange_rad, 'linear', 'extrap');
    
    % Get extinction difference properly
    measAbsCurve = (1-exp(-tempMeasEx*nanocrystalNumber(sizeIndex)/1));
    
    diffAbsCurve = 1 - (1 - measAbsCurve)./(1 - tempAbsCalc');
    
    diffExtCurve = -1/nanocrystalNumber(sizeIndex).*log(1-diffAbsCurve);
    
    dataToFitGroup{iSize} = diffExtCurve';
   
    frequencyGroup{iSize} = frequencyRange_rad;
end

%% Do for 2nd mode
%%% Subtraction approach didn't work so well, as fitting remainder
%%% influences main mode - need to do both in parallel

x0 = parametersToUse;

% Shift up resonant slope by measured ratios
x0(1) = x0(1)*nanmean(nanocrystal2ndFreqResonance_hz(sizesToUse)./nanocrystalFreqResonance_hz(sizesToUse));

lb = [x0(1)/2 0.001 1];
ub = [x0(1)*2 0.5 20];

f = @(x)(calculateerror_intensityfunction(x, dataToFitGroup, diameterDistGroup, ...
    numberGroup, frequencyGroup, 10^21, typeOfFunction));

options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 10^6, 'MaxIterations', 10^4);

solution = lsqnonlin(f, x0, lb, ub, options)

solution - x0
    
solution = [solution(1)*1000 solution(2:end-1)/10^6 solution(end)];

for iSize = 1:nSize
    subplot(2,nSize,iSize+nSize); hold on
    
    diametersToUse = diameterDistGroup{iSize};

    resonance_rad = 1./(diametersToUse/2) * solution(1);
    
    switch typeOfFit
        case '3Coeffs'
            intensityValues = solution(2);
            
        case 'IntensityPerSize'
            intensityValues = solution(2:end-1);
            
        case 'IntensityFunction'
            
            switch typeOfFunction
                case 'Linear'
                    % Diameter in nm
                    refSize = diametersToUse*10^9;

                case 'Quadratic'
                    % Area in nm^2
                    refSize =  4*pi*(diametersToUse/2).^2*10^18;

                case 'Cubic'
                    %volume in nm^3
                    refSize = 4/3*pi*(diametersToUse/2).^3*10^27;
            end
            
            intensityValues = refSize*solution(2);
    end
    
    % Check predicted cross section
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, resonance_rad, ...
        intensityValues, solution(end), numberGroup{iSize}, []);
    
    plot(frequencyRangePlot_rad/2/pi/10^9, (extinctionCrossSection)/10^-21, 'r--')
    
    sizeIndex = sizesToUse(iSize);
    
    plot(frequencyRangePlot_rad/2/pi/10^9, (extinctionCrossSection+mode1Data{iSize})/10^-21, 'r:')
    
    plot(frequencyGroup{iSize}/2/pi/10^9, dataToFitGroup{iSize}/10^-21, 'b')
     
    % Get full extinction curve
    measuredExtinctionCurve = interp1(curveFrequncy*2*pi, ExCrossSecCurve_m2{sizeIndex}, frequencyRangePlot_rad, 'linear', 'extrap');
    
    plot(frequencyRangePlot_rad/2/pi/10^9, measuredExtinctionCurve/10^-21, 'm-'); hold on
end
%%
% Calculates extinction curve difference from resonator. Multiple sizes, but all use equal coefficients
function sse = calculateerror_equalcoefficients(x, extinctionRef, diameters, number, ...
    frequencyRange, positionRef)
  
    % For solver
    % Expects all frequency values in radians and diameters in m
    
%     if length(x) ~= 2 + numberOfModes
%        error('wrong number of test parameters') 
%     end
    
    if iscell(extinctionRef)
       error('Not configured for cells') 
    end
    
    % Unpack test values
    slopeInds = find(positionRef == 1);
    slopeMultInds = find(positionRef == 2);
    chargeInds = find(positionRef == 3);
    QFInds = find(positionRef == 4); 
    
    testSlope = x(slopeInds)*1000;
    
    testSlopeMult = x(slopeMultInds);
    
    testCharge = x(chargeInds)/10^6;
    
    testQFs = x(QFInds);
    
    numberOfModes = length(QFInds);
    
    resonance_rad = 1./(diameters/2) * testSlope;
    resonanceBase_rad = resonance_rad;
    
    abosrbtion = zeros(length(frequencyRange),1);
    
    % Do for each mode and combine absorbtions
    for iMode = 1:numberOfModes
        if iMode > 1
            resonance_rad = resonanceBase_rad*testSlopeMult(iMode-1);
        end
        
        tempAbs = calculateresonantormixtureabsorbtion(frequencyRange, resonance_rad, ...
            testCharge, testQFs(iMode), number, []);
        
        abosrbtion = 1 - (1 - abosrbtion) .* (1 - tempAbs);
    end

    extinctionCrossSection = -1/sum(number).*log(1-abosrbtion);
    
    % Bottom for lsqnonlin - requires array of differences, then matches fmincon
    sse = (extinctionCrossSection-extinctionRef)/max(extinctionRef);
end

% Creates extinction curve from resonator. Multiple sizes and varying intensity - x(2)
function sse = calculateerror_varyingintensity(x, extinctionRef, diameters, number, ...
    frequencyRange, outPutScalar)

    % For solver
    % Expects, x(1): slope, x(2:end-1): intensity for each diameter (q/mass), x(end): gamma
    % Expects all frequency values in radians and diameters in m
    
    if length(x) ~= length(diameters)+2
       error('wrong number of test parameters') 
    end
    
    if iscell(extinctionRef)
       error('Not configured for cells') 
    end
    
    x = [x(1)*1000 x(2:end-1)/10^6 x(end)];
    
    resonance_rad = 1./(diameters/2) * x(1);
    
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRange, resonance_rad, ...
        x(2:end-1), x(end), number, []);

    if outPutScalar
        % Top for fmincon - requires single value
        sse = sum((conversion*(extinctionCrossSection-extinctionRef)).^2);
    else
        % For lsqnonlin - requires array of differences
        sse = (extinctionCrossSection-extinctionRef)/max(extinctionRef);
    end
end

% Creates extinction curve from resonator. Multiple sizes and intensity function
function sse = calculateerror_intensityfunction(x, extinctionRef, diameters, number, ...
    frequencyRange, typeOfFunction, positionRef)

    % For solver
    % Expects all frequency values in radians and diameters in m
    
%     if length(x) ~= 3
%        error('wrong number of test parameters') 
%     end
    
    % Decide if function is run with multiple sizes or not
    if iscell(extinctionRef)
       nSize = length(extinctionRef);
    else
       nSize = 1;
    end

    
    % Unpack test values
    slopeInds = find(positionRef == 1);
    slopeMultInds = find(positionRef == 2);
    chargeInds = find(positionRef == 3);
    QFInds = find(positionRef == 4); 
    
    testSlope = x(slopeInds)*1000;
    
    testSlopeMult = x(slopeMultInds);
    
    testCharge = x(chargeInds)/10^6;
    
    testQFs = x(QFInds);
    
    numberOfModes = length(QFInds);
    
    sse = [];
    
    for iSize = 1:nSize
        if nSize > 1
            freqToUse = frequencyRange{iSize};
            
            diameterToUse = diameters{iSize};
            
            numberToUse = number{iSize};
            
            extinctionRefToUse = extinctionRef{iSize};

        else
            freqToUse = frequencyRange;
            
            diameterToUse = diameters;
            
            numberToUse = number;
            
            extinctionRefToUse = extinctionRef;
        end
        
        switch typeOfFunction
            case 'Linear'
                % Diameter in nm
                refSize = diameterToUse*10^9;

            case 'Quadratic'
                % Area in nm^2
                refSize =  4*pi*(diameterToUse/2).^2*10^18;

            case 'Cubic'
                %volume in nm^3
                refSize = 4/3*pi*(diameterToUse/2).^3*10^27;
        end

        intensityValues = refSize * testCharge;

        resonance_rad = 1./(diameterToUse/2) * testSlope;
        resonanceBase_rad = resonance_rad;
        
        abosrbtion = zeros(length(frequencyRange),1);
    
        % Do for each mode and combine absorbtions
        for iMode = 1:numberOfModes
            if iMode > 1
                resonance_rad = resonanceBase_rad*testSlopeMult(iMode-1);
            end

            tempAbs = calculateresonantormixtureabsorbtion(freqToUse, resonance_rad, ...
                intensityValues, testQFs(iMode), numberToUse, []);

            abosrbtion = 1 - (1 - abosrbtion) .* (1 - tempAbs);
        end

        extinctionCrossSection = -1/sum(numberToUse).*log(1-abosrbtion);
        
%         [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(freqToUse, resonance_rad, ...
%             intensityValues, x(end), numberToUse, []);
        
        % Just concatanting SSE
        % Can be different lengths not sure if that is a problem...

        % normalize weighting, helps fit lower mag curves
        sse = [sse, (extinctionCrossSection-extinctionRefToUse)'/max(extinctionRefToUse)];
    end
    
    sse = sse';
end
