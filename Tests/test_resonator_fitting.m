clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Start with sizes measured on 50 GHz Bandwidth
sizesToUse = [1 4 6];    

nSize = length(sizesToUse);

simWithInterpSizeDist = 1; % If not set, uses original distribution for simulation

varyCharge = 0; % in simulation, using quadratic model?

sizeSteps = 1/10^9; % in m

fitRealData = 1; % if not set, fits simulated data

% fitWithSizeType = 'Single';
fitWithSizeType = 'Original';
% fitWithSizeType = 'FullInterp';

% typeOfFit = '3Coeffs';

% typeOfFit = 'IntensityPerSize';
    constrainToIncreasing = 0; %fmincon will be used with linear inequality

typeOfFit = 'IntensityFunction';
    typeOfFunction = 'Linear';
%     typeOfFunction = 'Quadratic';
%     typeOfFunction = 'Cubic';

useHandFreqLimits = 1;

if ~useHandFreqLimits
    frequencyRange_rad = (100:5:400)*10^9*2*pi;
end

frequencyRangePlot_rad = (50:1:450)*10^9*2*pi;

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
        frequencyRange_rad = (lowerFreqToFit(sizeIndex):5:upperFreqToFit(sizeIndex))*10^9*2*pi;
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

    % Get slope of main size
    resonantSlope = resonance_rad/(1/(nanocrystalSize_m(sizeIndex)/2));
    
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
    
    % Set up fit
    switch typeOfFit
    
        % Use solver for 3 coefficients
        case '3Coeffs'
            % Parameters rescalled to similar ranges
            x0 = [resonantSlope/1000, bulkCharge/sqrt(reducedMass)*10^6, bulkQF];
            
            lb = [resonantSlope/1000/2 0 0];
            ub = [resonantSlope/1000*2 20 20];
            
            f = @(x)(calculateerror_equalcoefficients(x, dataToFit, diametersToUse, ...
                numberToUse, frequencyRange_rad, 10^21));
    
        % Use solver allowing for varying intensity    
        case 'IntensityPerSize'
        %%% Soln values quite sensitive to initial parameter choices
            % Acceptable space deffinetely worth considering...
        x0 = [resonantSlope/1000, bulkCharge/sqrt(reducedMass)*10^6*ones(1,length(diametersToUse)),...
            bulkQF];
    
        % Test providing 2nd step of answer...
%         x0 = [resonantSlope/1000, 0.14:(1.3-0.14)/(length(diametersToUse)-1):1.3,...
%             bulkQF];
        
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
            numberToUse, frequencyRange_rad, 10^21, constrainToIncreasing));
    
        case 'IntensityFunction'
            
            x0 = [resonantSlope/1000, 0.01, bulkQF];

            lb = [resonantSlope/1000/2 0.001 1];
            ub = [resonantSlope/1000*2 0.5 20];
            
            f = @(x)(calculateerror_intensityfunction(x, dataToFit, diametersToUse, ...
                numberToUse, frequencyRange_rad, 10^21, typeOfFunction));
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
    
    solution = [solution(1)*1000 solution(2:end-1)/10^6 solution(end)];
    
    % Calculate parameters solved for
    resonance_rad = 1./(diametersToUse/2) * solution(1);
    
    if length(solution) == 3
        parametersToUseGroup(iSize, :) = solution;
    end
    
    % Intensity depends on type of solution
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
        intensityValues, solution(end), numberToUse, []);
    
    plot(frequencyRangePlot_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'r--')
end

% Run on combined!
x0 = mean(parametersToUseGroup);

lb = [resonantSlope/1000/2 0.001 1];
ub = [resonantSlope/1000*2 0.5 20];

f = @(x)(calculateerror_intensityfunction(x, dataToFitGroup, diameterDistGroup, ...
    numberGroup, frequencyGroup, 10^21, typeOfFunction));

options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 10^6, 'MaxIterations', 10^4);

solution = lsqnonlin(f, x0, lb, ub, options)

solution - x0
    
solution = [solution(1)*1000 solution(2:end-1)/10^6 solution(end)];

figure;

%Plot for each size
for iSize = 1:nSize
    subplot(1,nSize,iSize);
    
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
    
    plot(frequencyRangePlot_rad/2/pi/10^9, extinctionCrossSection/10^-21, 'r--')
    
end
%%
% Calculates extinction curve difference from resonator. Multiple sizes, but all use equal coefficients
function sse = calculateerror_equalcoefficients(x, extinctionRef, diameters, number, ...
    frequencyRange, conversion)
  
    % For solver
    % Expects, x(1): slope, x(2): intensity (q/mass), x(3): gamma
    % Expects all frequency values in radians and diameters in m
    
    if length(x) ~= 3
       error('wrong number of test parameters') 
    end
    
    if iscell(extinctionRef)
       error('Not configured for cells') 
    end
    
    x = [x(1)*1000 x(2)/10^6 x(3)];
    
    resonance_rad = 1./(diameters/2) * x(1);
    
    [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(frequencyRange, resonance_rad, ...
        x(2), x(3), number, []);

    % Apply conversion to error to keep it above solver tolerance.

    % Bottom for lsqnonlin - requires array of differences, then matches fmincon
    sse = conversion*(extinctionCrossSection-extinctionRef);
end

% Creates extinction curve from resonator. Multiple sizes and varying intensity - x(2)
function sse = calculateerror_varyingintensity(x, extinctionRef, diameters, number, ...
    frequencyRange, conversion, outPutScalar)

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
        sse = conversion*(extinctionCrossSection-extinctionRef);
    end
end

% Creates extinction curve from resonator. Multiple sizes and intensity function
function sse = calculateerror_intensityfunction(x, extinctionRef, diameters, number, ...
    frequencyRange, conversion, typeOfFunction)

    % For solver
    % Expects, x(1): slope, x(2): intensity (q/mass), x(3): gamma
    % Expects all frequency values in radians and diameters in m
    
    if length(x) ~= 3
       error('wrong number of test parameters') 
    end
    
    % Decide if function is run with multiple sizes or not
    if iscell(extinctionRef)
       nSize = length(extinctionRef);
    else
       nSize = 1;
    end
    
    x = [x(1)*1000 x(2)/10^6 x(end)];
    
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

        intensityValues = refSize * x(2);

        resonance_rad = 1./(diameterToUse/2) * x(1);

        [~, extinctionCrossSection] = calculateresonantormixtureabsorbtion(freqToUse, resonance_rad, ...
            intensityValues, x(end), numberToUse, []);
        
        % Just adding up SSE
        %%% may cause errors as each is a different lenght... could weight?
        
        % Bottom for lsqnonlin - requires array of differences
        sse = [sse, conversion*(extinctionCrossSection-extinctionRefToUse)'];
    end
    
    sse = sse';
end
