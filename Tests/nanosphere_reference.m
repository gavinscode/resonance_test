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

% Using text values for all but 8.4 (typo in text)
nanocrystalFreqBandwidth_hz = [31 45.70 34 40 17 40]*10^9;%+-

% Using text values for all
nanocrystalThetaEx_m2 = [0.7*10^-22 1.1*10^-21 1.6*10^-21 2.0*10^-21 4.1*10^-21 5.7*10^-21];

%Currently have lowest, centre and highest bin centres
%%% Get other bin centres an counts
nanocrystalSizeDist_binCentres_8 = [5.05 6.34 7.54 8.83 10.16 11.46]*10^-9;
nanocrystalSizeDist_Counts_8 = [2 11 20 9 6 2];

nanocrystalSizeDist_binCentres_10p4 = [5.44 7.26 9.12 10.98 12.70 14.56]*10^-9;
nanocrystalSizeDist_Counts_10p4 = [2 9 18 16 5 5];

nanocrystalSizeDist_binCentres_13 = [6.66 10.25 13.85 17.29 20.72]*10^-9;
nanocrystalSizeDist_Counts_13 = [6 17 49 12 1];

% Materials - CdSe is core, CdTe is shell
% From Savoit's app
% just using this (no CdTe fraction) fit's in middle of experimental points
% CdSeVelocity_mps = [3700, 1540];  % m/s 

% From Sun's paper 
% just using this (no, CdTe fraction, as Sun did), still towards higher points
% using both, still tends towards higher points
%CdSeVelocity_mps = [3570 1720];  % m/s 

%CdSeDensity_kgpm3 = 5810;        % kg/m^3  

%CdTeVelocity_mps = [3411, 1756]; % m/s

% From Adachi's book
% using both for [111] tends towards lower points
% using both for [100] fit's in middle
% Speeds for 100
CdSeVelocity_mps = [3430 1980];  % m/s 

CdSeDensity_kgpm3 = 5664;        % kg/m^3

CdTeVelocity_mps = [3020 1860];  % m/s 

CdTeDensity_kgpm3 = 5856;        % kg/m^3 


% From Calc in get_aperture file
apertureArea = 2.10*10^-4;

% General constants
VACCUM_PERMITIVITY = 8.854187817*10^-12; %C^2/(N.M^2)

LIGHT_SPEED = 299792458; % m/s