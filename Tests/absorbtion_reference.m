% First test, Liu ... Sun 2008 CdSe/CdTe type-II nanocrystals

useCompleteCurves = 1;

% Extinction cross sections curves stored in cell array
ExCrossSecCurve_m2 = cell(6,1);

if useCompleteCurves
    interpSpacing = 5;
    curveFrequncy = (100:interpSpacing:400)*10^9;
    
    cd('/Users/gavintaylor/Documents/Matlab/Git_versioned_April_27/Resonance_test_git/Data')
    
    % Use curves digitised from Sun Combined 2008 w/ WbPlotDigitizer
    data = readmatrix('Intensity_in.csv');
    % Extrapoalte each to have matched frequncy range
    I_source_mv = interp1(data(:,1)*10^9, data(:,2), curveFrequncy, 'linear', 'extrap');
    
    data = readmatrix('Intensity_13.csv');
    % Extrapoalte each to have matched frequncy range
    I_trans_13_mv = interp1(data(:,1)*10^9, data(:,2), curveFrequncy, 'linear', 'extrap');
    
    data = readmatrix('Ext_8.csv');
    ExCrossSecCurve_m2{1} = interp1(data(:,1)*10^9, data(:,2)*10^-21, curveFrequncy, 'linear', 'extrap');
    
    data = readmatrix('Ext_10p1.csv');
    ExCrossSecCurve_m2{3} = interp1(data(:,1)*10^9, data(:,2)*10^-21, curveFrequncy, 'linear', 'extrap');
    
    data = readmatrix('Ext_10p4.csv');
    ExCrossSecCurve_m2{4} = interp1(data(:,1)*10^9, data(:,2)*10^-21, curveFrequncy, 'linear', 'extrap');
    
    data = readmatrix('Ext_13.csv');
    ExCrossSecCurve_m2{6} = interp1(data(:,1)*10^9, data(:,2)*10^-21, curveFrequncy, 'linear', 'extrap');
    
    % Set frequency bounds for fitting - 1st mode
    %%% All multiuples of 5
    lowerFreqToFit_1st = [150, NaN NaN, 150, NaN, 110];
    
    upperFreqToFit_1st = [240, NaN NaN, 230, NaN, 210];
    
    lowerFreqToFit_2nd = [NaN, NaN NaN, 275, NaN, 200];
    
    upperFreqToFit_2nd = [NaN, NaN NaN, 400, NaN, 310];
else
    % Source intensity
    I_source_mv = [2.92, 2.58, 2.60, 2.33, 1.76, 1.18, 0.69, 0.42, 0.32, 0.34 0.34, 0.29, 0.22]; %mV

    % Transmitted intensity 13 nm - only have curve for 6th
    I_trans_13_mv = [2.87, 2.03, 1.72, 1.53, 1.19, 0.83, 0.49, 0.33, 0.25, 0.28, 0.27, 0.19, 0.16]; %mV

    %%% Added zeros to prevent fitted curve lifting at low values
    warning('zeros preprended below 100');

    % Frequncies measuresed
    curveFrequncy = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400]*10^9; % Hz

    % For 8 nm
    ExCrossSecCurve_m2{1} = [0, 0, 0, 0.12, 0.33, 0.32, 0.47, 0.67, 0.54, 0.14, 0.26, 0.15, 0.18 0.30, 0.47, 0.36]*10^-21; % m^2

    % For 10.1 nm
    ExCrossSecCurve_m2{3} = [0, 0, 0, 0.21, 0.37, 1.08, 1.62, 1.04, 0.25, 0.24, 0.90, 1.36, 3.34, 2.39, 1.52, 1.08]*10^-21; % m^2

    % For 10.4 nm
    ExCrossSecCurve_m2{4} = [0, 0, 0, 0.12, 0.38, 0.61, 0.94, 1.57, 1.66, 1.57, 1.65, 1.70, 1.97, 2.14, 2.18 1.98]*10^-21; % m^2

    % For 13 nm
    ExCrossSecCurve_m2{6} = [0, 0, 0, 0.57, 3.44, 5.49, 5.68, 5.1, 4.62, 4.49, 2.96, 2.42, 2.78, 3.41, 4.88, 3.38]*10^-21; % m^2
end