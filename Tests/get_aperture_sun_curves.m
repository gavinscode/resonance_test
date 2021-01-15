% From curves in SPIE paper
% Needed to determine apertue area
% Also wanted to check most plausible number of spheres for each size 

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
areaToUse = find((I_source_mv-I_trans_13_mv)./I_source_mv > 0.25);

% Just avg between 'reliable' points
avgArea = mean(apertureArea(areaToUse))

apertureRadius = sqrt(avgArea/pi)*1000



I_trans_10p4_calc = 10.^(ExCrossSecCurve_m2{sizesToUse(2)}*-numberSpheres(2)/avgArea) .* I_source_mv;

I_trans_8_calc = 10.^(ExCrossSecCurve_m2{sizesToUse(1)}*-numberSpheres(1)/avgArea) .* I_source_mv;



figure; subplot(1,2,1); hold on;

plot(curveFrequncy/10^9, I_source_mv, 'k')

plot(curveFrequncy/10^9, I_trans_13_mv, 'r')

plot(curveFrequncy(areaToUse)/10^9, I_trans_13_mv(areaToUse), 'rx')

plot(curveFrequncy/10^9, I_trans_10p4_calc, 'xm')

plot(curveFrequncy/10^9, I_trans_8_calc, 'xb')

title('Transmission intensity')


subplot(1,2,2); hold on;

plot(curveFrequncy/10^9, apertureArea, 'r')

plot(curveFrequncy(areaToUse)/10^9, apertureArea(areaToUse), 'rx')

line([curveFrequncy(1) curveFrequncy(end)]/10^9, [avgArea avgArea])

title('Aperture area')