% From curves in SPIE paper
% Needed to determine apertue area
% Also wanted to check most plausible number of spheres for each size 

% Above 300 GHz intensity still matches but extinction isn't so good.

clear; close all; clc

% Load data
nanosphere_reference

absorbtion_reference

%%% At resonance, should be a D^4 relation for cross-section
% - so 10^(-D^4*n/A)*I_source - did not check, but seems most likely for given order of spheres...

sizesToUse = [1 4 6];

% Smallest to largest - seems order is correct
numberSpheres = nanocrystalNumber(sizesToUse);
%numberSpheres = fliplr(numberSpheres);

apertureArea = ExCrossSecCurve_m2{sizesToUse(3)}*-numberSpheres(3)./log(I_trans_13_mv./I_source_mv);

%areaToUse = find(curveFrequncy/10^9 > 120 & curveFrequncy/10^9 < 230);
areaToUse = find((I_source_mv-I_trans_13_mv)./I_source_mv > 0.25 & I_source_mv > 0.5);

% Just avg between 'reliable' points
avgArea = mean(apertureArea(areaToUse))

apertureRadius = sqrt(avgArea/pi)*1000

%Calc intensity for other sizes
I_trans_10p4_calc = exp(ExCrossSecCurve_m2{sizesToUse(2)}*-numberSpheres(2)/avgArea) .* I_source_mv;

I_trans_8_calc = exp(ExCrossSecCurve_m2{sizesToUse(1)}*-numberSpheres(1)/avgArea) .* I_source_mv;


% Get intensity and exCross for 13 given avg area to test
I_trans_13_calc = exp(ExCrossSecCurve_m2{sizesToUse(3)}*-numberSpheres(3)/avgArea) .* I_source_mv;

exCrossSecCurve_13_calc = -avgArea/numberSpheres(3).*...
        log(I_trans_13_mv./I_source_mv);

    

figure; subplot(1,3,1); hold on;

plot(curveFrequncy/10^9, I_source_mv, 'k')

plot(curveFrequncy/10^9, I_trans_13_mv, 'r')

plot(curveFrequncy/10^9, I_trans_13_calc, 'rx')

plot(curveFrequncy(areaToUse)/10^9, I_trans_13_mv(areaToUse), 'ro')

plot(curveFrequncy/10^9, I_trans_10p4_calc, '-m')

plot(curveFrequncy/10^9, I_trans_8_calc, '-b')

title('Transmission intensity')

legend({'Source', '13 nm meas.', '13nm calculated', '13 nm used for avg. area', '10.4 calc', '8nm calc'})
ylabel('Intensity')


subplot(1,3,2); hold on;

plot(curveFrequncy/10^9, apertureArea, 'r')

plot(curveFrequncy(areaToUse)/10^9, apertureArea(areaToUse), 'ro')

line([curveFrequncy(1) curveFrequncy(end)]/10^9, [avgArea avgArea])

title('Aperture area (mv)')

legend({'Calcualted', 'Used in avg.', 'Avg value'})

ylabel('Aoertyre area (mm2)')


subplot(1,3,3); hold on;

plot(curveFrequncy/10^9, ExCrossSecCurve_m2{sizesToUse(3)}/10^-21, 'r')

plot(curveFrequncy/10^9, exCrossSecCurve_13_calc/10^-21, 'rx')

plot(curveFrequncy/10^9, ExCrossSecCurve_m2{sizesToUse(2)}/10^-21, '-m')

plot(curveFrequncy/10^9, ExCrossSecCurve_m2{sizesToUse(1)}/10^-21, '-b')

title('Extinction cross section')
legend({'13 nm meas' '13 nm calc', '10.4 meas', '8 nm meas'})
ylabel('Cross section (m2 * 10-21)')