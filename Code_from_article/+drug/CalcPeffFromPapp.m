function [P_eff] = CalcPeffFromPapp(MW, P_app)
% Calculate Peff from Papp according to Sjögren et al. 2013
% Input: Papp [cm/s] e.g. 7.6*10-6 cm/s, par

% 1) Isolate the membrane permeability (P_m) by removing the permeability by
%  the aqueous boundary layer (P_ABL)

% 1a) Calculate the diffusion coefficient in water
visc_factor  = 1;                               % Viscosity in ELF relative to saline
[D_w]        = drug.StokesEinstein(MW, visc_factor); % [dm2/h]

% Unit transformation from dm2/h to m2/s
D_w = D_w/3600/100; % 1/3600 --> dm2/h to dm2/s; 1/100 --> dm2/s to m2/s

% 1b) Define the thickness of ABL (Sjögren et al. 2013)
L_caco = 300e-6;    % 300 um
% L_caco = L_caco*10; % Transformation from m to dm

% 1c) Calculate P_ABL (Sjögren et al. 2013)
P_ABL = D_w/L_caco;  % m/s
P_ABL = P_ABL*100;   % [cm/s]

% 1d) Calculate P_m
P_m = 1/((1/P_app) - (1/P_ABL)); % [cm/s]

% 2) Transformation to Peff using empirical relaionship from Sjögren et al.
%    2013. 
%    NOTE! Unit is 10-4 cm/s; i.e. 5e-6 cm/s should enter the eq as 0.05
P_m_eq = P_m*1e4;
log_p_eff = 0.9762*log10(P_m_eq) + 1.4012; % [10-4 cm/s]
P_eff = (10^log_p_eff)*1e-4; % back to 10-6 cm/s

end