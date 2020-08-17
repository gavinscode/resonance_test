% For radius
%%% All data for spherical, could consider variation from spherical as well.

% A: From Parupudi 2017 - H3N2, possibly inactivated
diameterMean_range_m_A = [114 153]*10^-9; % Table 1 - means

diameterRange_m_A = [70 300]*10^-9; % Figure 1 - extremes from all

%%% Re take std as percentage of mean given range %%%
    % Probably more consistent
% C: From Ruigrok 1984 for X49 - H3N2 x H1N1
diameterMean_m_C = 140*10^-9;

diameterSD_m_C = 12*10^-9;

% D: From Ruigrok 1985 - For X31 (his other papers are in these ranges)
diameterMean_range_m_D = [123 138]*10^-9;

diameterSD_range_m_D = [12 43]*10^-9;

% E: From Moules 2011
diameterMean_range_m_E = [88.2 - 99.0]*10^-9;

diameterSD_range_m_E = [8.5 - 18.0]*10^-9;

% For density - convert from MDa to N
massMean_range_kg_A = [265 306]*1.6605*10^-21; % From table 2

massSD_range_kg_A = [66 147]*1.6605*10^-21; % From table 2
% A - Supp lists monomer value as 206 [144 268], unclear what this represents?

massMean_kg_C = 161*1.6605*10^-21;

massSD_kg_C = 17*1.6605*10^-21;

massRange_kg_C = [90 240]*1.6605*10^-21;

%%% reduced mass fraction...?

% For sound velocity
% Assume ratio of 0.5 from longitudinal to transverse



% For charge
%%% Where are there opposing charges: core, capsid, lipid membrane?


% For fracture stress
%%% Need information for ruptering capsids

% B: From Li 2011 - for lipid membrane
fractureForceRange_N_B = [0.25 2000]*10^-12; %N

% Force/Area given 30 nm radius AFM tip
fractureStress_Pa_B = fractureForceRange/(pi*(30*10^-9)^2); %Pa
