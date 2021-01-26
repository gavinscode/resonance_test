clear; close all; clc

% Load data
nanosphere_reference

% Interpolate Core Size
sizeWithDist = [1 4 6];

scaleCoreSize = 1;

coresToGet = find(isnan(nanocrystalCore_m));

nanocrystalCore_m(coresToGet) = interp1(nanocrystalSize_m(sizeWithDist), nanocrystalCore_m(sizeWithDist), ...
    nanocrystalSize_m(coresToGet), 'linear');

% Calc q and Q values from original to get rough idea on relationship
sizesToUse = 1:6;    

qEstimate = zeros(length(sizesToUse),1);

QEstimate = zeros(length(sizesToUse),1);

% Get sphere parameters
sphereArea = 4*pi*(nanocrystalSize_m(sizesToUse)/2).^2;

sphereVolume = 4/3*pi*(nanocrystalSize_m(sizesToUse)/2).^3;

sizeSteps = 0.1/10^9;

avgSphereRadius = zeros(length(sizeWithDist),1);

avgSphereArea = zeros(length(sizeWithDist),1);
 
avgSphereVolume = zeros(length(sizeWithDist),1);

qEstimateVolDist = zeros(length(sizeWithDist),1);

for iSize = 1:length(sizesToUse)
    sizeIndex = sizesToUse(iSize);
   
    QEstimate(iSize) = nanocrystalFreqResonance_hz(sizeIndex)/nanocrystalFreqBandwidth_hz(sizeIndex);
    
    resonance_rad = nanocrystalFreqResonance_hz(sizeIndex)*2*pi;

    if any(sizeWithDist == sizeIndex)
        sizeInd = find(sizeWithDist == sizeIndex);
        
        % Get size dist
        diameterDist = nanocrystalSizeDist{sizeIndex, 1};

        countDist = nanocrystalSizeDist{sizeIndex, 2};

        % Coarse and fine distribution means are the same
        [countDist, diameterDist] = interpolatescaleddistribution(countDist, diameterDist, sizeSteps);

        avgSphereRadius(sizeInd) = wmean(diameterDist, countDist);

        avgSphereArea(sizeInd) = wmean(4*pi*(diameterDist/2).^2, countDist);

        avgSphereVolume(sizeInd) = wmean(4/3*pi*(diameterDist/2).^3, countDist);
        
        % Also calc q given this volume distribution
        reducedMass = zeros(length(diameterDist),1);
        
        for jDiameter = 1:length(diameterDist)
            if scaleCoreSize
                % Scale core diameter given distribution
                %%% Note, assumse core directly scales with total diamter, may not be true
                coreDiameterScaled_m = nanocrystalCore_m(sizeIndex)*diameterDist(jDiameter)/ ...
                    nanocrystalSize_m(sizeIndex);
            else
                coreDiameterScaled_m = nanocrystalCore_m(sizeIndex);
            end

            coreVolume = 4/3*pi*(coreDiameterScaled_m/2).^3;

            totalVolume = 4/3*pi*(diameterDist(jDiameter)/2).^3;

            coreMass = coreVolume*CdSeDensity_kgpm3;

            % Take difference in volume for shell mass
            if totalVolume - coreVolume > 0
                shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;
            else
               % If core larger than shell, don't let it go negative
               shellMass = 0;
            end
        
            reducedMass(jDiameter) = coreMass*shellMass/(coreMass + shellMass);
        end
        
        countDist(reducedMass == 0) = [];
        
        reducedMass(reducedMass == 0) = [];
        
        qEstimateVolDist(sizeInd) = wmean(sqrt(nanocrystalThetaEx_m2(sizeIndex) .* resonance_rad .* ...
            reducedMass.*VACCUM_PERMITIVITY.*LIGHT_SPEED./QEstimate(iSize)), countDist');
    end
    
    coreVolume = 4/3*pi*(nanocrystalCore_m(sizeIndex)/2).^3;

    totalVolume = 4/3*pi*(nanocrystalSize_m(sizeIndex)/2).^3;

    coreMass = coreVolume*CdSeDensity_kgpm3;

    % Take difference in volume
    shellMass = (totalVolume - coreVolume)*CdTeDensity_kgpm3;

    reducedMass = coreMass*shellMass/(coreMass + shellMass);

%     systemSpring = resonance_rad^2*reducedMass;
% 
%     systemDamp = resonance_rad*reducedMass/QEstimate(iSize);
    
    qEstimate(iSize) = sqrt(nanocrystalThetaEx_m2(sizeIndex) * resonance_rad * ...
        reducedMass*VACCUM_PERMITIVITY*LIGHT_SPEED./QEstimate(iSize)); 
end

% Test size comparisons
figure;
subplot(1,3,1); hold on
plot(nanocrystalSize_m*10^9, nanocrystalSize_m*10^9, 'b-x'); 
plot(nanocrystalSize_m(sizeWithDist)*10^9, avgSphereRadius*10^9, 'rx'); 
title('Diameter')

subplot(1,3,2); hold on
plot(sphereArea*10^18, sphereArea*10^18, 'b-x'); 
plot(sphereArea(sizeWithDist)*10^18, avgSphereArea*10^18, 'rx'); 
title('Area')

subplot(1,3,3); hold on
plot(sphereVolume*10^27, sphereVolume*10^27, 'b-x'); 
plot(sphereVolume(sizeWithDist)*10^27, avgSphereVolume*10^27, 'rx'); 
title('Volume')


figure;
subplot(3,2,1); hold on
plot(nanocrystalSize_m(sizesToUse)*10^9, QEstimate, 'bx-')
plot(nanocrystalSize_m(sizeWithDist)*10^9, QEstimate(sizeWithDist), 'bo')
plot(avgSphereRadius*10^9, QEstimate(sizeWithDist), 'rx-')
title('Quality factor');
ylim([0 15]); xlim([5 15]);

subplot(3,2,2); hold on
plot(nanocrystalSize_m(sizesToUse)*10^9, qEstimate/(1.602176634*10^-19), 'bx-')
plot(nanocrystalSize_m(sizeWithDist)*10^9, qEstimate(sizeWithDist)/(1.602176634*10^-19), 'bo')
plot(avgSphereRadius*10^9, qEstimateVolDist/(1.602176634*10^-19), 'rx-')
title('Charge (in e) - not fits are on blue');
ylim([0 500]); 
xlim([0 15]);

% Note this is e as function of nanometers...

% Use fit 2nd order polynomial
warning('Note - line only uses sizes with distributions')
opts = fitoptions('poly1', 'Lower', [-Inf 0], 'Upper', [Inf 0]);
fittedLine = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineH = plot(fittedLine,'k');

opts = fitoptions('poly1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedLine = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineHFree = plot(fittedLine,'k:');

opts = fitoptions('poly2', 'Lower', [-Inf 0 0], 'Upper', [Inf 0 0]);
fittedQuad = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'poly2', opts)
quadH = plot(fittedQuad,'b');

opts = fitoptions('exp1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedExp = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'exp1', opts)
expH = plot(fittedExp, 'g');

opts = fitoptions('power1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedPower = fit([nanocrystalSize_m']*10^9, [qEstimate]/(1.602176634*10^-19), 'power1', opts)
powerH = plot(fittedPower, 'r');

legend([lineH, quadH, expH, powerH], {'Linear', 'Quad', 'Exp', 'Power'})

subplot(3,2,4); hold on
plot(sphereArea*10^18, qEstimate/(1.602176634*10^-19), 'bx-')
plot(sphereArea(sizeWithDist)*10^18, qEstimate(sizeWithDist)/(1.602176634*10^-19), 'bo')
plot(avgSphereArea*10^18, qEstimateVolDist/(1.602176634*10^-19), 'rx-')

title('Charge (in e)');
ylim([0 500]); 
xlim([0 600]);

opts = fitoptions('poly1', 'Lower', [-Inf 0], 'Upper', [Inf 0]);
fittedLine = fit([sphereArea']*10^18, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineH = plot(fittedLine,'k');

opts = fitoptions('poly1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedLine = fit([sphereArea']*10^18, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineHFree = plot(fittedLine,'k:');

subplot(3,2,6); hold on
plot(sphereVolume*10^27, qEstimate/(1.602176634*10^-19), 'bx-')
plot(sphereVolume(sizeWithDist)*10^27, qEstimate(sizeWithDist)/(1.602176634*10^-19), 'bo')
plot(avgSphereVolume*10^27, qEstimateVolDist/(1.602176634*10^-19), 'rx-')

title('Charge (in e)');
ylim([0 500]); 
xlim([0 1500]);

opts = fitoptions('poly1', 'Lower', [-Inf 0], 'Upper', [Inf 0]);
fittedLine = fit([sphereVolume']*10^27, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineH = plot(fittedLine,'k');

opts = fitoptions('poly1', 'Lower', [-Inf -Inf], 'Upper', [Inf Inf]);
fittedLine = fit([sphereVolume']*10^27, [qEstimate]/(1.602176634*10^-19), 'poly1', opts)
lineHFree = plot(fittedLine,'k:');