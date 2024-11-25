function [Cd,vg,ne,pa,Dmol] = ParametersA(da,po)
% po=unit density; da=aerodynamic diameter
% Cd=cummingham slip correction factor; vg=gravitational settling velocity
% of a particle [m/s]; Dmol=particle diffusion constant [m2/s]; 
% ne=viscosity of air [kg/m/s]; pa=density of air (kg/m3)
% Use SI-units consistently

da = da.*10^-6; % Transformation from um to m

k  = 1.38064852*10^-23 ; % boltzmann constant [m^2·kg/(s^2·K)] 
T  = 273.15+37.5;        % absolute temp, 37.5 degree Celcius (range 37-37 degree C)
ne = 1.9224364E-5;       % viscosity of air at 37.5 dgr C [kg/m/s] (1.846*10^-5 kg/m/s at 300K) 
g  = 9.81;               % gravitational acceleration [m/s^2]
pa = 1.1372;             % density of air 37.5 degr C [kg/m3]; http://www.gribble.org/cycling/air_density.html
lamda=0.066*10^-6;       % mean free path of air molecule; 0.066 um, http://myweb.uiowa.edu/tpeters/IH1/Aerosols/AerosolFormulae.pdf [m]

% Cd from: Lee et al. and https://en.wikipedia.org/wiki/Cunningham_correction_factor

Cd   = 1 + (lamda./da).*(2.514+0.8*exp(-0.55*(da./lamda))); % Cummingham slip correction factor
vg   = (po.*((da.^2).*g.*Cd)./(18*ne));                     % Gravitational settling velocity of a particle
Dmol = (k*T.*Cd)./(3*pi*ne.*da);                            % Brownian diffusion coefficient  [m2/s]

end
