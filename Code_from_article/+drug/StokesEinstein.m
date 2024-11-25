function [D_ES] = StokesEinstein(MW, visc_factor)
% Calculation of diffusion coefficient [dm2/h] using Stokes-Einstein
% equation.
% Input needed: MW [g/mol or Da] and visc_factor; viscosity of ELF 
% compared to water. Example: 
% visc_factor = 1 --> equal viscosity
% visc_factor = 10 --> 10-fold higher viscosity in ELF

% Einstein-Stokes requires radius of molecule (r), which can be estimated
% from molecular weight (MW) according to:
% r = MW*p(1)+p(2); % Radius in nm

% Parameter values obtained from optimization
p = [0.000834373134971   0.178704501830042];

% Einstein-Stokes -note that the viscosity is 7-10-fold higher in ELF relative 
% to saline (healthy volunteers) according to ref b.

T = 273.15 + 37; % Absolute temp in Kelvin
vis = visc_factor*0.6913*1e-3; % 0.6913 cPa, cP = 1 mPa新 = 0.001 Pa新 = 0.001 N新搶?2 = 0.001 kg搶?1新?1.
kB = 1.38064852*1e-23; %m2 kg s-2 K-1
r = MW*p(1)+p(2); % Radius in nm
r = r*1e-9; 

D_ES = kB*T/(6*pi*vis*r); % [m2/s]

% Unit transformation from m2/s to dm2/h

D_ES = D_ES*3600*100;

end

% References:
% a) Plant Cell Biology: From Astronomy to Zoology; 2009; Randy Wayne. Chapter 3, p. 57
% b) http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3114535/

% Reference values of r: 0.51 and 0.76 nm, approximate radius for MW = 400 and 700 Da, resp. (ref a) 