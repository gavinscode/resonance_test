%% To determine q (charge distribution) values

% From Fig. 3B of Yang 2016 
    % (y) absorbtion measured in imageJ for each Ghz
    % (x) frequncy measured for peak and zero
    % Peak is actually at 8.3??? Changed from 8.22
calcFrequencies_hz = [6 7 8 8.3 9 10 11 12 13.25 14]*10^9;

calcFrequencies_rad = calcFrequencies_hz*2*pi;

calcAbsorbtions = [9.7 16.4 21.1 21.5 20.1 13.8 6.4 1.7 0.1 0.5]/100;

if length(calcFrequencies_hz) ~= length(calcAbsorbtions)
   error('Incorectly copied values') 
end

peakInd = 4;

% from Ellison et al. 1996 pg 240 Pottel 1980, 25oC
permitivityValues = [73.11, 73.19, 73.13, 72.71, 72.63, 72.25, 71.92, 71.50, ... 
    71.21, 70.34, 69.45, 69.74, 69.24, 68.78, 68.69, 67.86, 67.13, 66.77, ... 
    65.38, 63.76, 63.04, 61.74, 62.04, 61.10, 59.46, 58.49, 58.18, 55.78, ... 
    55.27, 54.03, 53.47, 52.63]; 

permitivityFrequenceis = [5.306, 5.323, 5.433, 5.536, 5.638, 5.853, 6.000, ...
    6.145, 6.300, 6.729, 6.850, 6.958, 7.267, 7.406, 7.681, 7.850, 8.243, ...
    8.579, 8.979, 9.516, 10.010, 10.150, 10.230, 10.450, 11.320, 11.730, ...
    12.000, 12.770, 13.140, 13.380, 13.820, 14.230];

if length(permitivityValues) ~= length(permitivityFrequenceis)
   error('Incorectly copied values') 
end

calcRelativePermitivity = interp1(permitivityFrequenceis*10^9, permitivityValues, ...
    calcFrequencies_hz, 'linear','extrap');

% Calcualte values for q
% Eqn 7: q omitted as it is solved for and E is ommited as it will cancel later
calcAmplitudeWOq = 1./(reducedMass_kg*...
        sqrt((measuredResonance_rad^2 - calcFrequencies_rad.^2).^2 + ...
        (measuredResonance_rad*calcFrequencies_rad/measuredQualityFactor).^2));
 
% Eqn 9 - now omitting q^2   
calcAvgAbsorbtion = measuredResonance_rad.*calcFrequencies_rad.^2*...
    reducedMass_kg.*calcAmplitudeWOq.^2/(2*measuredQualityFactor);

% E^2 is ommited as will cancel later
calcPowerFlux = 0.5*sqrt(calcRelativePermitivity)*VACCUM_PERMITIVITY*LIGHT_SPEED;

% Eqn 10, q^2 is omitted as it will be solved for
calcThetaAbs1 = calcAvgAbsorbtion./calcPowerFlux;

% Eqn 13
calcThetaAbs2 = -1/(virionDensity*channelLength_m).*log(1-calcAbsorbtions);

% Calculate Q from equvielent between Theta Abs
calcQ = sqrt(calcThetaAbs2./calcThetaAbs1);

% Relative to orignal charge distrubiton
calcQ/providedChargeDistribution;

% Compare to eqn 14 - simplified by assuming freq at resonance
sqrt(calcThetaAbs2(peakInd)*measuredResonance_rad*reducedMass_kg*...
    sqrt(calcRelativePermitivity(peakInd))*VACCUM_PERMITIVITY*LIGHT_SPEED...
    /measuredQualityFactor)/providedChargeDistribution;

% Just using q from peak reduces quality and increases bandwidth a bit...
%calcQ(:) = providedChargeDistribution;

% Interpolate for charge distribution up to top freq calculated
qInterpolated = interp1(calcFrequencies_hz, calcQ, testFrequncies_hz, ...
    'linear','extrap');

% From Murrey et al 2006 Eqn 23, sets practical limit for Q 
qLimit = sqrt(0.15*calcFrequencies_rad(peakInd)*calcThetaAbs2(peakInd)*sqrt(calcRelativePermitivity(peakInd))*...
    VACCUM_PERMITIVITY*LIGHT_SPEED*reducedMass_kg);

centreQinCharge = calcQ(peakInd)/(1.602176634*10^-19)

providedChargeDistribution = calcQ(peakInd);

limitInCharge = qLimit/(1.602176634*10^-19)

% Quick plot to test results
figure; hold on; 

plot(testFrequncies_hz/10^9, qInterpolated/providedChargeDistribution)

plot(calcFrequencies_hz/10^9, calcQ/providedChargeDistribution, 'ro')

plot(calcFrequencies_hz([1 end])/10^9, [1 1]*qLimit/providedChargeDistribution)

% Interpolate permitivity
relativePermitivtyInterpolated = interp1(permitivityFrequenceis*10^9, permitivityValues, ...
    testFrequncies_hz, 'linear','extrap');