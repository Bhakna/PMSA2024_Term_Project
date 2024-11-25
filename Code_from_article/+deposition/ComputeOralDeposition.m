function [E_0] = ComputeOralDeposition(Q, d_um, Dmol)
% This function calculates oral deposition as described by Cheng
% et al 2003 "Aerosol Deposition in the Extrathoracic Region"
% Input: d_um [um], diffusion coefficient Dmol [m2/s], 
% inspiratory/expiratory flow rate Q [mL/min] 

% Original paper by Cheng et  al 1996 uses the following units
% D  [cm2/s]
% da [um]
% Q  [L/min]

Q_eq = Q*60/1000; % unit transformation from mL/s to L/min

% The diffusion coefficient is transformed from m2/s to cm2/s inside of the
% formulae

E_0 = 1 - exp(-0.000278*Q_eq*d_um.^2 - 20.4*(10000*Dmol).^0.66*Q_eq^-0.31); % Oral breathing

end

