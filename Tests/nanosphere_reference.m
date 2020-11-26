% First test, Liu ... Sun 2008 CdSe/CdTe type-II nanocrystals

%%% Also aim to get data from Murray and other nanosphere L1 papers

nanocrystalSize_nm = [8 8.4 10.1 10.4 11.6 13]; %diameter?

nanocrystalNumber = [3.6*10^15 NaN NaN 1.4*10^16 NaN 1.4*10^16];

nanocrystalCore_nm = [4.3 NaN NaN 4.3 NaN 5,3]; %diameter?

nanocrystalFreqResolution_GHz = [50 10 10 50 10 50];

%Using plot measurments for 8, 8.4 (missing in text), 10.4
nanocrystalFreqResonance_GHz = [198 241 175 221 170 165];

% Did not measure from plots, 8.4 is missing in text
nanocrystalFreqBandwidth_GHz = [31 NaN 34 40 17 40]%+-

%Using plot measurments for 8, 10.4, 
nanocrystalThetaEx_m2 = [6.74*10^-22 1.1*10^-21 1.6*10^-21 1.67*10^-21 4.1*10^-21 5.7*10^-21];

%Currently have lowest, centre and highest bin centres
%%% Get other bin centres an counts
nanocrystalSizeDist_binCentres_8 = [(4.47+5.52)/2 NaN (7.04+8.13)/2 NaN NaN (14.04+15.41)/2];
nanocrystalSizeDist_Counts_8 = [NaN NaN NaN NaN NaN NaN];

nanocrystalSizeDist_binCentres_10p4 = [(4.75+6.2)/2 NaN (8.39+9.93)/2 NaN NaN (14.04+15.41)/2];
nanocrystalSizeDist_Counts_10p4 = [NaN NaN NaN NaN NaN NaN];

nanocrystalSizeDist_binCentres_13 = [(5.41+8.08)/2 NaN (12.31+15.29)/2 NaN (19.43+22.04)/2];
nanocrystalSizeDist_Counts_13 = [NaN NaN NaN NaN NaN];

CdSeVl = 3570; % m/s - core material?

CdSeVt = 1720; % m/s - core material?