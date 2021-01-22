% First test, Liu ... Sun 2008 CdSe/CdTe type-II nanocrystals

%%% Also get data from Murray and other nanosphere L1 papers

nanocrystalSize_m = [8 8.4 10.1 10.4 11.6 13]*10^-9; 

% Weird order, but correct given aperture calcs
nanocrystalNumber = [3.6*10^15 NaN NaN 1.4*10^16 NaN 1.6*10^16];

%%% Core variation is indicated in methods paper
    % Assumed values for 8.4 10.1 and 11.6
nanocrystalCore_m = [4.3 NaN NaN 4.3 NaN 5.4]*10^-9;

% Note: 50 is resolution bandwidth/FWHM
    % 10 is spectral resolution, unsure of relation to gaussian
freqResolution_Ghz = [50 10 10 50 10 50];

% Using plot measurments for 8, 8.4 (missing in text), 10.4
% nanocrystalFreqResonance_hz = [198 241 175 221 170 165]*10^9;

% Using text values for all but 8.4 (typo in text)
nanocrystalFreqResonance_hz = [210 241 175 200 170 165]*10^9;

nanocrystal2ndFreqResonance_hz = [380  NaN 322 360 271 250]*10^9;

nanocrystal3rdFreqResonance_hz = [NaN NaN NaN NaN 333 370]*10^9;

% Using text values for all but 8.4 (typo in text)
nanocrystalFreqBandwidth_hz = [31 45.70 34 40 17 40]*10^9; % Note, this is HW

nanocrystal2ndFreqBandwidth_hz = [50 NaN 32 70 15 NaN]*10^9;

nanocrystal3rdFreqBandwidth_hz = [NaN NaN NaN NaN 15 60]*10^9;

% Double as HW provided in text - FW required for getting systemQ
nanocrystalFreqBandwidth_hz = nanocrystalFreqBandwidth_hz*2;

% Using text values for all
nanocrystalThetaEx_m2 = [0.7 1.1 1.6 2.0 4.1 5.7]*10^-21;

% Size ditributions stored in cell array - only have for 1st, 4th and 6th
% [diameter, counts]
nanocrystalSizeDist = cell(6,2);

nanocrystalSizeDist{1,1} = [5.05 6.34 7.54 8.83 10.16 11.46]*10^-9; %bin center diameter, m
nanocrystalSizeDist{1,2} = [2 11 20 9 6 2]; %count

nanocrystalSizeDist{4,1} = [5.44 7.26 9.12 10.98 12.70 14.56]*10^-9;
nanocrystalSizeDist{4,2} = [2 9 18 16 5 5];

nanocrystalSizeDist{6,1} = [6.66 10.25 13.85 17.29 20.72]*10^-9;
nanocrystalSizeDist{6,2} = [6 17 49 12 1];

% Materials - CdSe is core, CdTe is shell
% From Savoit's app
% just using this (no CdTe fraction) fit's in middle of experimental points
%CdSeVelocity_mps = [3700, 1540];  % m/s %ratio: 2.4

% From Sun's paper 
% just using CdSe (no, CdTe fraction, as Sun did), tends towards higher points
% using both, maintains similar fit
CdSeDensity_kgpm3 = 5810;        % kg/m^3  

%CdSeVelocity_mps = [3570 1720];  % m/s %ratio:2.07

%CdTeVelocity_mps = [3411, 1756]; % m/s %ratio: 1.9

% Calculated using spherical average of stiffness constants in Adachi 2005
CdSeVelocity_mps = [3625 1733];  % m/s %ratio: 2.1

CdTeVelocity_mps = [3280 1623];  % m/s %ratio: 2.0

% W limits from Murray site - for testing
% Ratio < 2 tends to shift right, ratio > 2 tends to shift left
%CdSeVelocity_mps = [3501 1493]; % low low %ratio: 2.3 %shift left
%CdSeVelocity_mps = [3501 1825]; % low high %ratio: 1.9 % shift right
%CdSeVelocity_mps = [3750 1493]; % high low %ratio: 2.51 % shift left
%CdSeVelocity_mps = [3750 1825]; % high high %ratio: 2.05 % shift right
%CdTeVelocity_mps = CdSeVelocity_mps;

CdTeDensity_kgpm3 = 5856;        % kg/m^3 

% From Calc in get_aperture file
apertureArea = 2.11*10^-4;

% General constants
VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

LIGHT_SPEED = 299792458; % m/s