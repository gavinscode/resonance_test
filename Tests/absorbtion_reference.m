% First test, Liu ... Sun 2008 CdSe/CdTe type-II nanocrystals

% Source intensity
I_source_mv = [2.92, 2.58, 2.60, 2.33, 1.76, 1.18, 0.69, 0.42, 0.32]; %mV

% Transmitted intensity 13 nm - have curve for 1 size
I_trans_13_mv = [2.87, 2.03, 1.72, 1.53, 1.19, 0.83, 0.49, 0.33, 0.25]; %mV

% Frequncies measuresed
curveFrequncy = [100, 125, 150, 175, 200, 225, 250, 275, 300]*10^9; % Hz

% Extinction cross sections curves stored in cell array
ExCrossSecCurve_m2 = cell(6,1);

% For 8 nm
ExCrossSecCurve_m2{1} = [0.12, NaN, 0.32, NaN, 0.67, NaN, 0.14, NaN, 0.15]*10^-21; % m^2

% For 10.4 nm
ExCrossSecCurve_m2{4} = [0.12, NaN, 0.61, NaN, 1.57, NaN, 1.57, NaN, 1.70]*10^-21; % m^2

ExCrossSecCurve_m2{6} = [0.51, 3.36, 5.43, 5.64, 5.08, 4.59, 4.43, 2.99, 2.40]*10^-21; % m^2
