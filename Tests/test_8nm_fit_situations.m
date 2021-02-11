clear; clc; close all

% Load data
nanosphere_reference

absorbtion_reference

% Only working with 8m for now
sizesToUse = 1;

nSize = length(sizesToUse);

frequencyRangePlot_rad = (100:1:500)*10^9*2*pi;

fitFreqStep = 1;

% From combined paper fig. 2f, average top and bottom
    % for m^2 as fn of freq in radians
    % background from 10.4 is way to high for 8nm
backgroundEstimate = (20*10^-22/(0.6*10^12*2*pi)^2 + 110*10^-22/(1.35*10^12*2*pi)^2)/2;

% For two 'pure' curves
if 0
    frequenciesToFit_rad = [106:fitFreqStep:128 258:fitFreqStep:276]*10^9*2*pi;
    
    guessFreqs = [119 267]*2*pi*10^9;

    guessQF = [5 15];

    % Really not sure how to get intensity from curve, for sphere it's charge/sqrt(mass)
    guessIntensities = [1 1]*8.5*10^-7;
    
    numberOfModes = 1;
    
    useBackground = 0;
    
    useConstraints = 0;
    
    % When these are fitted with others, they are lifted substantially by
    % other modes
end

% For four first order curves + 2nd and 3rd modes
if 1
    frequenciesToFit_rad = [100:fitFreqStep:500]*10^9*2*pi;
    
    guessFreqs = [119 165 198 267]*2*pi*10^9;

    guessQF = ones(1,length(guessFreqs))*5;

    guessIntensities = ones(1,length(guessFreqs))*8.5*10^-7;
    
    numberOfModes = 5;
    
    guessQFScale = ones(1, numberOfModes-1);
    
    guessVl_Vt = (CdSeVelocity_mps + CdTeVelocity_mps)/2;
    
    useBackground = 0;
    
    useConstraints = 0;
end

% For six first order curves
    % Note this is just to test, expecting two to be 2nd order
if 0
    frequenciesToFit_rad = [100:fitFreqStep:314]*10^9*2*pi;
    
    guessFreqs = [119 165 198 228 267 301]*2*pi*10^9;

    guessQF = ones(1,length(guessFreqs))*5;

    guessIntensities = ones(1,length(guessFreqs))*8.5*10^-7;
    
    numberOfModes = 1;
    
    useBackground = 0;
    
    useConstraints = 0;
end

% For many first order curves
    % Note this is just to test
    % Fits fairly well, but needs to use constraints to keep everything in place
if 0
    frequenciesToFit_rad = [100:fitFreqStep:500]*10^9*2*pi;
    
    % May be something at 385 or could just be broad peak - works better w/o
    % Needs point above 500, doesn't compensate with background
    guessFreqs = [119 165 198 228 267 296 334 370 431 460 500]*2*pi*10^9; 

    guessQF = ones(1,length(guessFreqs))*5;

    guessIntensities = ones(1,length(guessFreqs))*8.5*10^-7;
    
    numberOfModes = 1;
    
    useBackground = 0;
    
    useConstraints = 1;
        liftUpLast = 1;
end

nResonators = length(guessFreqs);

for iSize = 1:nSize
    sizeIndex = sizesToUse(iSize);

    % Get curve
    measuredExtinctionCurve = interp1(ExCrossSecCurve_m2{sizeIndex,1}*2*pi, ExCrossSecCurve_m2{sizeIndex,2},...
        frequencyRangePlot_rad, 'linear', 'extrap');
    
    measuredExtinctionCurve(measuredExtinctionCurve < 0) = 0;
    
    % Take points to fit
    [~, fitInds] = intersect(frequencyRangePlot_rad, frequenciesToFit_rad);
    
    % Rescaled
    x0 = [guessFreqs/10^11 guessIntensities*10^6 guessQF];
    
    if useBackground
        x0 = [x0 backgroundEstimate*10^46/1000];
    end
    
    if numberOfModes > 1
       x0 = [x0 guessQFScale guessVl_Vt]; 
    end
    
    if useConstraints 
        lb = [guessFreqs*0.9/10^11, 0.01*ones(1,nResonators), ones(1,nResonators)];

        % Lift up for last
        if liftUpLast
            ub = [[guessFreqs(1:end-1)*1.1 600*10^9*2*pi]/10^11, 2*ones(1,nResonators), 50*ones(1,nResonators)];
        else
           ub = [guessFreqs*1.1/10^11, 2*ones(1,nResonators), 50*ones(1,nResonators)]; 
        end
        
        if useBackground
            lb = [lb 0];
            
            ub = [ub backgroundEstimate*10^46];
        end
        
        if numberOfModes > 1
           error('No constraints set yet') 
        end
    else
       lb = [];
       
       ub = [];
    end
        
    f = @(x)(calculateerror_seperate(x, measuredExtinctionCurve(fitInds), frequenciesToFit_rad, ...
        nResonators, useBackground, numberOfModes));

    options = optimoptions('lsqnonlin','Display','iter', 'FunctionTolerance', 1e-10, ...
            'MaxFunctionEvaluations', 10^6, 'MaxIterations', 20);
        
    solution = lsqnonlin(f, x0, lb, ub, options)
    
    solution - x0
    
    %%% Surely a function that can be made to do the unpacking and curve calculations,
        %%% it's the same as in error calculation function
    
    solvedResonance = solution(1:nResonators)*10^11;
    
    solvedIntensity = solution(nResonators+1:2*nResonators)/10^6;
    
    solvedQF = solution(2*nResonators+1:3*nResonators);
    
    [solvedResonance'/2/pi/10^9, solvedQF' solvedIntensity'*10^6]'
    
    if useBackground
        
        solvedBackground = solution(3*nResonators+1)/10^46;
        
        solvedBackground/backgroundEstimate
        
        indexOffset = 1;
    else
        indexOffset = 0;
    end
    
    if numberOfModes > 1
        solvedQFscale = [1 solution(3*nResonators+indexOffset+1:3*nResonators+indexOffset+numberOfModes-1)]
        
        solvedSoundSpeeds = solution(3*nResonators+indexOffset+numberOfModes:3*nResonators+indexOffset+numberOfModes+1)
        
        % Get frequencies then ratio to first
        freqs = calcualtesphereresonance(10/10^9/2, 'sph', 1, numberOfModes-1, ...
            solvedSoundSpeeds(1), solvedSoundSpeeds(2), 5*10^9, 10^6, 0)*2*pi;
        
        solvedFreqRatios = freqs./freqs(1)
    else
        solvedQFscale = 1;
        
        solvedFreqRatios = 1;
    end
    
    %%% Modify for extra modes, as in original
    % Simulate output
    numberOfSpheres = ones(nResonators, 1);
    
    extinctionCrossSection = zeros(length(frequencyRangePlot_rad),1);
    
    % Note that solution input is rescaled
    for iMode = 1:numberOfModes
        [~, tempEx] = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, solvedResonance.*solvedFreqRatios(iMode), ...
                solvedIntensity, solvedQF.*solvedQFscale(iMode), numberOfSpheres, 1, []);
            
        extinctionCrossSection = extinctionCrossSection + tempEx;    
    end
    %%% Should get function to correct for this itself..    
    extinctionCrossSection = extinctionCrossSection*nResonators;    
        
    if useBackground
        background = solvedBackground*frequencyRangePlot_rad.^2;
        
        extinctionCrossSection = extinctionCrossSection + background';
    end
    
    % plot main curve
    figure;  
    subplot(2,1,1); hold on
    plot(frequencyRangePlot_rad/2/pi/10^9, measuredExtinctionCurve*10^21, 'm')
    
    plot(frequenciesToFit_rad/2/pi/10^9, measuredExtinctionCurve(fitInds)*10^21, 'b.')
    
    plot(frequencyRangePlot_rad/2/pi/10^9, extinctionCrossSection*10^21, 'b')
    
    resColors = jet(nResonators);
    
    for jResonator = 1:nResonators
        individualCrossSection = zeros(length(frequencyRangePlot_rad),1);
        
        for iMode = 1:numberOfModes
            [~, tempEx] = calculateresonantormixtureabsorbtion(frequencyRangePlot_rad, solvedResonance(jResonator)*solvedFreqRatios(iMode), ...
                solvedIntensity(jResonator), solvedQF(jResonator)*solvedQFscale(iMode), 1, 1, []);
            
            individualCrossSection = individualCrossSection + tempEx;
        end
        
        plot(frequencyRangePlot_rad/2/pi/10^9, individualCrossSection*10^21, 'color', resColors(jResonator, :))
    end
    
    if useBackground
        plot(frequencyRangePlot_rad/2/pi/10^9, background*10^21, 'k')
    end
    
    % Plot residiual
    subplot(2,1,2)
    
    plot(frequenciesToFit_rad/2/pi/10^9, (measuredExtinctionCurve(fitInds)-extinctionCrossSection(fitInds)')*10^21, 'k')
    xlim([100 500])
end

% Gets extinction curve error for resonator collection. Treats each
% resonator seperately and only first mode
function sse = calculateerror_seperate(x, extinctionRef, ...
    frequencyRange, nResonators, includeBackground, nModes)

    if iscell(extinctionRef)
       error('Not configured for cells') 
    end
    
    % Number doesn't mean much for just extinctiong,
        % but needs to be large to prevent rounding errors
    numberOfSpheres = ones(nResonators, 1);

    % Note that x input is rescaled
    testResonance = x(1:nResonators)*10^11;
    
    testIntensity = x(nResonators+1:2*nResonators)/10^6;
    
    testQF = x(2*nResonators+1:3*nResonators);
    
    if includeBackground
        testBackground = x(3*nResonators+1)/10^46;
        
        indexOffset = 1;
    else
        indexOffset = 0;
    end
    
    if nModes > 1
        QFscale = [1 x(3*nResonators+indexOffset+1:3*nResonators+indexOffset+nModes-1)];
        
        soundSpeeds = x(3*nResonators+indexOffset+nModes:3*nResonators+indexOffset+nModes+1);
        
        % Get frequencies then ratio to first
        % Would be better to precalculate scattered interpolant.
        freqs = calcualtesphereresonance(10/10^9/2, 'sph', 1, nModes-1, ...
            soundSpeeds(1), soundSpeeds(2), 5*10^9, 10^6, 0)*2*pi;
        
        freqRatios = freqs./freqs(1);
    else
        freqRatios = 1;
        
        QFscale = 1;
    end
    
    extinctionCrossSection = zeros(length(frequencyRange),1);
    
    for iMode = 1:nModes
        [~, tempEx] = calculateresonantormixtureabsorbtion(frequencyRange, testResonance.*freqRatios(iMode), ...
            testIntensity, testQF*QFscale(iMode), numberOfSpheres, 1, []);
        
        extinctionCrossSection = extinctionCrossSection + tempEx;
    end

    % Scale up, as function calcualtes weighted average of all resonators
    extinctionCrossSection = extinctionCrossSection*nResonators;
    
    if includeBackground
        background = testBackground*frequencyRange.^2;
        
        extinctionCrossSection = extinctionCrossSection + background';
    end
    
    % For lsqnonlin - requires array of differences
    sse = (extinctionCrossSection-extinctionRef')/max(extinctionRef);
end
