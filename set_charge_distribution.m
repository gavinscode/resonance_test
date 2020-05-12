%% To determine q (charge distribution) values

%%% 10.5 may be a bit off as bend in curve
% From Fig. 3B of Yang 2016
calcFrequencies_hz = [6 8.22 10.5 11.25 13]*10^9;

calcFrequencies_rad = calcFrequencies_hz*2*pi;

calcAbsorbtions = [0.1 0.21 0.1 0.05 0.001];

% from Ellison et al. 1996 pg 240 Pottel 1980, 25oC
% Exact freqs used are: 
% 6.000 8.243 10.450 11.320 13.140
calcRelativePermitivity = [71.92 67.13 61.10 59.46 55.27]; % Relative, 

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
sqrt(calcThetaAbs2(2)*measuredResonance_rad*reducedMass_kg*...
    sqrt(calcRelativePermitivity(2))*VACCUM_PERMITIVITY*LIGHT_SPEED...
    /measuredQualityFactor)/providedChargeDistribution;

% Just using q from peak reduces quality and increases bandwidth a bit...
%calcQ(:) = providedChargeDistribution;

% Interpolate for charge distribution up to top freq calculated
qInterpolated = interp1(calcFrequencies_hz, calcQ, testFrequncies_hz, ...
    'linear','extrap');

% From Murrey et al 2006 Eqn 23, sets practical limit for Q 
qLimit = sqrt(0.15*calcFrequencies_rad(2)*calcThetaAbs2(2)*sqrt(calcRelativePermitivity(2))*...
    VACCUM_PERMITIVITY*LIGHT_SPEED*reducedMass_kg);

centreQinCharge = calcQ(2)/(1.602176634*10^-19)

limitInCharge = qLimit/(1.602176634*10^-19)

% Quick plot to test results
figure; hold on; 

plot(testFrequncies_hz/10^9, qInterpolated/providedChargeDistribution)

plot(calcFrequencies_hz/10^9, calcQ/providedChargeDistribution, 'ro')

plot(calcFrequencies_hz([1 end])/10^9, [1 1]*qLimit/providedChargeDistribution)

% Interpolate permitivity
relativePermitivtyInterpolated = interp1(calcFrequencies_hz, calcRelativePermitivity, ...
    testFrequncies_hz, 'linear','extrap');

%% Messy test for IR - but x17,000 gain!!! - probably not realistic
%%% Note, this does not really need to be done, will be single q value for particle

IR_WL = 1671*10^-9; % Middle 3rd water window
IR_freq_Rad = LIGHT_SPEED/IR_WL*2*pi;
% Assume absorbtion is simliar, so ThetaAbs2 will match as well
IR_permit = 1.77;

% This ends up super small - seems fair
IR_Amp = 1./(reducedMass_kg*...
        sqrt((measuredResonance_rad^2 - IR_freq_Rad.^2).^2 + ...
        (measuredResonance_rad*IR_freq_Rad/measuredQualityFactor).^2));

% This ends up relatively small - seems fair
IR_AvgAbs = measuredResonance_rad.*IR_freq_Rad.^2*...
    reducedMass_kg.*IR_Amp.^2/(2*measuredQualityFactor);

% This ends up a bit smaller - expected due to permitivity
IR_Flux = 0.5*sqrt(IR_permit)*VACCUM_PERMITIVITY*LIGHT_SPEED;

% This ends up relatively small - seems fair
IR_ThetaAbs1 = IR_AvgAbs/IR_Flux;

% Segelstein gives similar absorbtion coefficient for water between third water window and 8 Ghz, 
% Does this mean absorbtion ratio is similiar? If so, is cross-section equal?
IR_ThetaAbs2 = 2.5*10^-13;

% Because thetaAbs1 decrease but thetaAbs2 stays constant, q increases a lot...
IR_calcQ = sqrt(IR_ThetaAbs2./IR_ThetaAbs1);
    
IR_calcQ/providedChargeDistribution;