% For radius and mass (convert from MDa to kg)
%%% All data for spherical, could consider variation from spherical as well.

% A: From Parupudi 2017 - H3N2, possibly inactivated
diameterMean_m_A = 114*10^-9; % Table 1 - mean cryoTEM
    % Distribution skewed low has top tail
    
diameterPeak_m_A = 120*10^-9; % Figure S2 - measured from plot
    % Bin size is ~20, so overlaps mean
    
diameterRange_m_A = [60 160]*10^-9; % Figure S2 - measured from plot

diameterTopTail_m_A = 200*10^-9; % Figure S2 - measured from plot

%%% Compare density to general
massPeak_range_kg_A = [229 249]*1.6605*10^-21; % From table 2, mean values probably effected by aggregation

massCV_range_A = [66/265 147/306]; % From table 2, SD on Mean.

massNominal_kg_A = 206*1.6605*10^-21; % From supp.

massNominal_range_kg_A = [144 269]*1.6605*10^-21;

%%% This results in a very high density

% B: From Ruigrok 1984 for X49 - H3N2 x H1N1
diameterMean_m_B = 140*10^-9; 

diameterCV_B = 12*10^-9/diameterMean_m_B;

massMean_kg_B = 161*1.6605*10^-21;
    % Distribution skewed low has top tail

massPeak_kg_B = 155*1.6605*10^-21;
    % Bin size is 10, so nearly at mean
    
massCV_B = 17*1.6605*10^-21/massMean_kg_B;

massRange_kg_B = [90 210]*1.6605*10^-21;

massTopTail_kg_B = 240*1.6605*10^-21;

%%% This results in a very high density

% C: From Ruigrok 1985 - For X31 (his other papers are in these ranges)
diameterMean_range_m_C = [124 138]*10^-9;

diameterCV_range_C = [24/130 43/135 42/127 36/138 23/129 14/127 12/124 17/129 18/131 23/129]; % From table 1

diameterCV_range_C = [min(diameterCV_range_C) max(diameterCV_range_C)];

% D: From Moules 2011
diameterMean_range_m_D = [88.2 99.0]*10^-9;

diameterCV_range_D = [8.9/88.2 18.0/99.0 8.5/91.7 12.2/96.0]; % From supp.

diameterCV_range_D = [min(diameterCV_range_D) max(diameterCV_range_D)];

% E: From Harris 2006
diameterMean_m_E = 120*10^-9; 

diameterRange_m_E = [84 170]*10^-9;

% Core mass fraction - Schlze 1973
% Note, this is of RNP, not including the protein capsid
coreFraction_range = [0.8 1.1]*10/100;

% For sound velocity
% Assume ratio of 0.5 from longitudinal to transverse



% For charge
%%% Where are there opposing charges: core, capsid, lipid membrane?


% For fracture stress
%%% Need information for ruptering capsids

% From Li 2011 - for lipid membrane
fractureForceRange_N = [0.25 2000]*10^-12; %N

% Force/Area given 30 nm radius AFM tip
fractureStress_Pa = fractureForceRange/(pi*(30*10^-9)^2); %Pa
