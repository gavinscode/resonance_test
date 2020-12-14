% First test, s toolLiu ... Sun 2008 CdSe/CdTe type-II nanocrystals

%%% Also get data from Murray and other nanosphere L1 papers

nanocrystalSize_m = [8 8.4 10.1 10.4 11.6 13]*10^-9; 

% Weird order, but correct given aperture calcs
nanocrystalNumber = [3.6*10^15 NaN NaN 1.4*10^16 NaN 1.4*10^16];

%%% Core variation is indicated in methods paper
nanocrystalCore_m = [4.3 NaN NaN 4.3 NaN 5.3]*10^-9;

freqResolution_hz = [50 10 10 50 10 50]*10^9;

%Using plot measurments for 8, 8.4 (missing in text), 10.4
%%% Updates seem correct in Fig. 2a
nanocrystalFreqResonance_hz = [198 241 175 221 170 165]*10^9;

% Did not measure from plots, 8.4 is missing in text
%%% Can be read from Fig. 2a to confirm... (clearer in comb. paper)
nanocrystalFreqBandwidth_hz = [31 NaN 34 40 17 40]*10^9;%+-

%Using plot measurments for 10.4, 
nanocrystalThetaEx_m2 = [0.7*10^-22 1.1*10^-21 1.6*10^-21 1.67*10^-21 4.1*10^-21 5.7*10^-21];

%Currently have lowest, centre and highest bin centres
%%% Get other bin centres an counts
nanocrystalSizeDist_binCentres_8 = [(4.47+5.52)/2 NaN (7.04+8.13)/2 NaN NaN (10.99+12.04)/2]*10^-9;
nanocrystalSizeDist_Counts_8 = [NaN NaN NaN NaN NaN NaN];

nanocrystalSizeDist_binCentres_10p4 = [(4.75+6.2)/2 NaN (8.39+9.93)/2 NaN NaN (14.04+15.41)/2]*10^-9;
nanocrystalSizeDist_Counts_10p4 = [NaN NaN NaN NaN NaN NaN];

nanocrystalSizeDist_binCentres_13 = [(5.41+8.08)/2 NaN (12.31+15.29)/2 NaN (19.53+22.04)/2]*10^-9;
nanocrystalSizeDist_Counts_13 = [NaN NaN NaN NaN NaN];

% From Savoit's tool - closer to experiment, should check...
CdSeVelocity_mps = [3700, 1540];  % m/s 
%From Sun's paper 
%CdSeVelocity_mps = [3570 1720];  % m/s 

CdSeDensity_kgpm3 = 5810;

CdTeVelocity_mps = [3411, 1756]; % m/s

CdTeDensity_kgpm3 = []; 