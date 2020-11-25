%% To determine q (charge distribution) values

% Given paramters measurment from paper, Q will vary. But this is a fudge
% to plot the initial absorbtion.

if length(measuredFrequencies_hz) ~= length(measuredAbsorbtions)
   error('Incorectly copied values') 
end

[~, peakInd] = max(measuredAbsorbtions);

if length(permitivityValues) ~= length(permitivityFrequenceis)
   error('Incorectly copied values') 
end

calcRelativePermitivity = interp1(permitivityFrequenceis*10^9, permitivityValues, ...
    measuredFrequencies_hz, 'linear','extrap');

% Calcualte values for q
% Eqn 7: q omitted as it is solved for and E is ommited as it will cancel later
calcAmplitudeWOq = 1./(reducedMass_kg*...
        sqrt((measuredResonance_rad^2 - measuredFrequencies_rad.^2).^2 + ...
        (measuredResonance_rad*measuredFrequencies_rad/measuredQualityFactor).^2));
 
% Eqn 9 - now omitting q^2   
calcAvgAbsorbtion = measuredResonance_rad.*measuredFrequencies_rad.^2*...
    reducedMass_kg.*calcAmplitudeWOq.^2/(2*measuredQualityFactor);

% E^2 is ommited as will cancel later
calcPowerFlux = 0.5*sqrt(calcRelativePermitivity)*VACCUM_PERMITIVITY*LIGHT_SPEED;

% Eqn 10, q^2 is omitted as it will be solved for
calcThetaAbs1 = calcAvgAbsorbtion./calcPowerFlux;

% Eqn 13
calcThetaAbs2 = -1/(virionDensity*channelLength_m).*log(1-measuredAbsorbtions);

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
qInterpolated = interp1(measuredFrequencies_hz, calcQ, testFrequncies_hz, ...
    'linear','extrap');

% From Murrey et al 2006 Eqn 23, sets practical limit for Q 
qLimit = sqrt(0.15*measuredFrequencies_rad(peakInd)*calcThetaAbs2(peakInd)*sqrt(calcRelativePermitivity(peakInd))*...
    VACCUM_PERMITIVITY*LIGHT_SPEED*reducedMass_kg);

centreQinCharge = calcQ(peakInd)/(1.602176634*10^-19)

providedChargeDistribution = calcQ(peakInd);

limitInCharge = qLimit/(1.602176634*10^-19)

% Quick plot to test results
figure; hold on; 

plot(testFrequncies_hz/10^9, qInterpolated/providedChargeDistribution)

plot(measuredFrequencies_hz/10^9, calcQ/providedChargeDistribution, 'ro')

plot(measuredFrequencies_hz([1 end])/10^9, [1 1]*qLimit/providedChargeDistribution)

% Interpolate permitivity
relativePermitivtyInterpolated = interp1(permitivityFrequenceis*10^9, permitivityValues, ...
    testFrequncies_hz, 'linear','extrap');