% From curves in SPIE paper

% Source values
I_source = [2.92, 2.58, 2.60, 2.33, 1.76, 1.18, 0.69, 0.42, 0.32]; %mV

% For 13 nm - transmitted
I_trans_13 = [2.87, 2.03, 1.72, 1.53, 1.19, 0.83, 0.49, 0.33, 0.25]; %mV

ExCrossSec_13 = [0.51, 3.36, 5.43, 5.64, 5.08, 4.59, 4.43, 2.99, 2.40]; %x10^-21 m^2

% For 10.4 nm
ExCrossSec_10p4 = [0.12, NaN, 0.61, NaN, 1.57, NaN, 1.57, NaN, 1.70]; %x10^-21 m^2

% For 8 nm
ExCrossSec_8 = [0.12, NaN, 0.32, NaN, 0.67, NaN, 0.14, NaN, 0.15]; %x10^-21 m^2

testFrequncy = [100, 125, 150, 175, 200, 225, 250, 275, 300]; % GHz

% Smallest to largest - seems order is correct
numberSpheres = [3.6*10^15, 1.4*10^16, 1.6*10^16];
%numberSpheres = fliplr(numberSpheres);

apertureArea = ExCrossSec_13*10^-21*-numberSpheres(3)/log(I_trans_13/I_source);

% Just avg between 'reliable' points
avgArea = mean(apertureArea(find(testFrequncy > 120 & testFrequncy < 230)))

apertureRadius = sqrt(avgArea/pi)*1000



I_trans_10p4_calc = 10.^(ExCrossSec_10p4*10^-21*-numberSpheres(2)/avgArea) .* I_source;

I_trans_8_calc = 10.^(ExCrossSec_8*10^-21*-numberSpheres(1)/avgArea) .* I_source;



figure; subplot(1,2,1); hold on;

plot(testFrequncy, I_source, 'k')

plot(testFrequncy, I_trans_13, 'r')

plot(testFrequncy, I_trans_10p4_calc, 'xm')

plot(testFrequncy, I_trans_8_calc, 'xb')

title('Transmission intensity')


subplot(1,2,2); hold on;

plot(testFrequncy, apertureArea, 'r')

line([testFrequncy(1) testFrequncy(end)], [avgArea avgArea])
title('Aperture area')