function [DEP_i, SED_i, DIF_i] = DepositionProbabilityBreathHold(eps_b, D_i, da, t_b, phi)
% Calculation of deposition probability for the tree during breath-hold
% eps_b is specifically calculated for breath-hold
% sedimentation equation for breath-hold is taken from Yeh and Schum 1980
% as the one by Thomas gives complex numbers for large t_b

D_i = D_i.*10^-2; % Transformation from cm to m

% Input parameters for the diffusion
da    = da.*10^-6; % Transformation from um to m

k     = 1.38064852*10^-23 ; % boltzmann constant [m^2·kg/(s^2·K)] 
T     = 273.15+37.5;        % absolute temp, 37.5 degree Celcius (range 37-37 degree C)
ne    = 1.9224364E-5;       % viscosity of air at 37.5 dgr C [kg/m/s] (1.846*10^-5 kg/m/s at 300K) 
g     = 9.81;               % gravitational acceleration [m/s^2]
pa    = 1.1372;             % density of air 37.5 degr C [kg/m3]; http://www.gribble.org/cycling/air_density.html
lamda = 0.066*10^-6;        % http://myweb.uiowa.edu/tpeters/IH1/Aerosols/AerosolFormulae.pdf


% Cd from: Lee et al. and https://en.wikipedia.org/wiki/Cunningham_correction_factor
Cd = 1 + (lamda./da).*(2.514+0.8*exp(-0.55*(da./lamda))); % Cummingham slip correction factor

% Calculations

[n, m] = size(eps_b); % n=number of gen; m=number of particle sizes

DEP_i = zeros(n,m);
SED_i = zeros(n,m);
DIF_i = zeros(n,m);

% Probability of impaction during breath-hold = 0

% Sedimentation equation for breath-hold is taken from Yeh and Schum 1980
for j=1:m
    for i=2:n

        r_i   = D_i(i)/2;  % radius of tube
        r_i_p = da(j)/2;   % radius of particle
        SED_i(i,j) = 1 - exp( (-4*g*Cd(j)*(r_i_p^2)*t_b*cos(phi(i)))/(9*pi*ne*r_i) );     


    end
end

% Diffusion during breath-hold from Schum and Yeh 1980 
for j = 1:m
    for i = 2:n
        r_i = D_i(i)/2;
        r_i_p = da(j)/2; % radius of particle
        DIF_i(i,j) = 1 - exp(-5.784*k*T*Cd(j)*t_b/(6*pi*ne*r_i_p*(r_i^2))); 
    end
end



