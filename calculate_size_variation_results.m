%close all;
clear;
clc

% Load and/or calculate values
set_parameters

set_charge_distribution

% Clear ones that are being set/tested
clear diameter_m

clear reducedMass_kg

clear resonanceFrequency_hz, clear resonantFrequency_rad

clear qualityFactor, clear systemSpring, clear systemDamp

%clear providedChargeDistribution, 
%% Set/test paramters specific to this analysis
% From Lamb and Krug in Fields Virology 2001
averageDiameter_virus_m = 100/10^9;
rangeDiameter_virus_m = [80 120]/10^9;

% From Ruigrok 1984 for X49 - H3N2 x H1N1
% Mass looks to have somewhat skewed normal distrubiotn, with tail above 200 nm
%%% Could be better to model as skew disturbiton
averageMass_X49_kg = 161*1.6605*10^-21; % Convert from MDa

stdMass_X49_kg = 17*1.6605*10^-21;

rangeMass_X49_kg = [90 240]*1.6605*10^-21; %240

% Check expected mass range to 3 STD
expectedRangeMass_X49_kg = [-3 3]*stdMass_X49_kg;

% Ratio bewteen actual and expected mass range 
    %- mass range is ~50% larger than expected
(rangeMass_X49_kg-averageMass_X49_kg)./expectedRangeMass_X49_kg;

averageDiameter_X49_m = 140/10^9;

stdDiameter_X49_m = 12/10^9;

% Calculate volume and density
averageVolume_X49_m3 = 4/3*pi*(averageDiameter_X49_m/2)^3;

% Note std of volume requires adding to mean, then subtracting.
stdVolume_X49_m3 = 4/3*pi*((averageDiameter_X49_m+stdDiameter_X49_m)/2)^3 - ...
    averageVolume_X49_m3;

virusDensity = averageMass_X49_kg/averageVolume_X49_m3;

% Std ratio between measured and calculated 
    %-std of mass is x~2.5 lower than expected from radius
%%% Check that scaling std in this way is correct    
stdMass_X49_kg/(stdVolume_X49_m3*virusDensity);

% Ratio between reported volume and 100 nm sphere 
    %-paper uses x~3 expectation for 100 nm virus!
4/3*pi*(averageDiameter_virus_m/2)^3*virusDensity/averageMass_X49_kg;

% Check std scales to fit 100 nm diameter range
expectedStdDiameter_virus_nm = stdDiameter_X49_m*averageDiameter_virus_m/averageDiameter_X49_m;

expectedRangeDiameter_virus_nm = [-3 3]*expectedStdDiameter_virus_nm;

% Ratio bewteen 'accepted' and expected mass range 
    %- Range is ~20% smaller than expected 
(rangeDiameter_virus_m-averageDiameter_virus_m)./expectedRangeDiameter_virus_nm;

% take diameter from 75 to 125 to catch tails - start w/5 nm steps
valuesDiameter_virus_m = (75:5:125)'/10^9;

% Get distrubiton of diameters from normal - funciton does not like nm
sizeFrequncy = normpdf(valuesDiameter_virus_m*10^9, averageDiameter_virus_m*10^9, ...
    expectedStdDiameter_virus_nm*10^9);

sizeFrequncy = sizeFrequncy/sum(sizeFrequncy);

warning('mass scaled back to that of paper')
massScaleUp = averageMass_X49_kg/(4/3*pi*(averageDiameter_virus_m/2)^3*virusDensity);

reducedMassScale = 0.09; %Unclear why reduced, core is taken to be 10% volume

valuesMass_virus_kg = 4/3*pi*(valuesDiameter_virus_m/2).^3*virusDensity*massScaleUp;

valuesReducedMass_virus_kg = valuesMass_virus_kg*reducedMassScale;

%%% May wish to try with orignal reduced mass for comparison

% From Savoit - Lamb calculator using spherical, L=1, N=0 (1 on page)
    % https://saviot.cnrs.fr/lamb/index.html
    % Set Vl=1434 and Vt=717 from trial and error so 100 nm freq = 8.2 GHz
    % Not much lower than lower end on 2009 paper
valuesFrequncy_virus_hz = [10.943; 10.259; 9.655; 9.119; 8.639; 8.207; 7.816; ...
    7.461; 7.137; 6.839; 6.566]*10^9;

valuesFrequncy_virus_rad = valuesFrequncy_virus_hz*2*pi;

if length(valuesFrequncy_virus_hz) ~= length(valuesDiameter_virus_m)
   error('Check resonance list') 
end

%[valuesDiameter_virus_m, valuesFrequncy_virus_hz]

QValuesToTest = (1:5)*2;
%% Simulate results at each q value

%%% Currently just using charge distribution (q) from original
 % providedChargeDistribution 

% Store results QxDxF matrix, as simulation done for each 
analyticalAmplitude = zeros(length(QValuesToTest), length(valuesDiameter_virus_m), length(testFrequncies_hz)); 
 
analyticalAbsorbtion = zeros(length(QValuesToTest), length(valuesDiameter_virus_m), length(testFrequncies_hz)); 

analyticalStress = zeros(length(QValuesToTest), length(valuesDiameter_virus_m), length(testFrequncies_hz)); 

%warning('just using provided charge distibution (q)')
%qToUse = providedChargeDistribution*ones(length(testFrequncies_hz),1);
qToUse = qInterpolated;

figure;

for iQVal = 1:length(QValuesToTest)
    
   for jDiam = 1:length(valuesDiameter_virus_m)
       
       % Set system parameters for given diameter
       systemSpring = valuesFrequncy_virus_rad(jDiam)^2*valuesReducedMass_virus_kg(jDiam);
       
       systemDamp = valuesFrequncy_virus_rad(jDiam)*valuesReducedMass_virus_kg(jDiam)/...
           QValuesToTest(iQVal);
       
       for kFreq = 1:length(testFrequncies_hz)
            % Eqn 7
            analyticalAmplitude(iQVal, jDiam, kFreq) = qToUse(kFreq)./(valuesReducedMass_virus_kg(jDiam)*...
                sqrt((valuesFrequncy_virus_rad(jDiam)^2 - testFrequncies_rad(kFreq)^2)^2 + ...
                (valuesFrequncy_virus_rad(jDiam)*testFrequncies_rad(kFreq)/QValuesToTest(iQVal)).^2));

            % Eqn 9
            analyticalPower = valuesFrequncy_virus_rad(jDiam)*...
                testFrequncies_rad(kFreq)^2 * valuesReducedMass_virus_kg(jDiam)*...
                analyticalAmplitude(iQVal, jDiam, kFreq)^2 / (2*QValuesToTest(iQVal));

            powerFlux = 0.5*sqrt(relativePermitivtyInterpolated(kFreq))*VACCUM_PERMITIVITY*LIGHT_SPEED;

            % Eqn 10
            theta_absorbtion = analyticalPower/powerFlux;

            % Eqn 13 - solution for absorbtion
            analyticalAbsorbtion(iQVal, jDiam, kFreq) = (1-exp(-theta_absorbtion*...
                virionDensity.*sizeFrequncy(jDiam)*channelLength_m));

            % Eqn 11
            analyticalStress(iQVal, jDiam, kFreq) = 2*systemSpring*analyticalAmplitude(iQVal, jDiam, kFreq)/...
                (0.58*pi*(valuesDiameter_virus_m(jDiam)/2)^2);
           
       end
       % Plot results for each diameter
       subplot(5,3,(iQVal-1)*3+1); hold on
       plot(testFrequncies_hz/10^9, permute(analyticalAmplitude(iQVal, jDiam, :), [3 2 1])*10^12);

       subplot(5,3,(iQVal-1)*3+2); hold on
       plot(testFrequncies_hz/10^9, permute(analyticalAbsorbtion(iQVal, jDiam, :), [3 2 1]));
       
       subplot(5,3,(iQVal-1)*3+3); hold on
       plot(testFrequncies_hz/10^9, permute(analyticalStress(iQVal, jDiam, :), [3 2 1]));
   end
   
   % Plot weighted sum - only makes sense for absorbtion
   subplot(5,3,(iQVal-1)*3+2); hold on
   temp = permute(analyticalAbsorbtion(iQVal, :, :), [2 3 1]);
   plot(testFrequncies_hz/10^9, sum(temp), 'k');
end