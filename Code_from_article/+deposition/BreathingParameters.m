function [breath] = BreathingParameters()
% A function for defining breathing parameters needed for calculating
% regional drug deposition in man

% Assign delay volume (V_D), tidal volume (V_T) and bolus volume
% (V_B) in mL
breath.V_D = 0;                   % Delay volume [mL] 
breath.V_T = 1000;                % Tidal volume [mL]  
breath.V_B = 1000;                % Bolus volume [mL] 

if (breath.V_D + breath.V_B) > breath.V_T
    error('The sum of the delay and bolus volume cannot exceed the tidal volume')
end

% Define breathing conditions
% --> define length of inspiration (t_in) and expiration (t_exp),
%     breathing paus (t_b), tidal volyme (V_T), breathing frequency (f_br)
    
breath.f_br   = 15;                    % Breathing frequency (1/min)  
breath.fr_in  = 0.5;                   % Fraction of breath as inspiratory
breath.t_b    = 0;                     % Breath-hold time [s]

% Define functional residual capacity (FRC)
breath.FRC = 3000;                     % Functional residual capacity, FRC [mL]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculations of t_in, t_exp, Q_in and Q_exp based on user-input 

breath.t_in   = breath.fr_in*(1/breath.f_br);        % Inspiratory time [min]
breath.t_exp  = (1-breath.fr_in)*(1/breath.f_br);    % Expiratory time  [min]

breath.Q_in   = breath.V_T/breath.t_in/60;           % Inspiratory flow rate [cm3/s]
breath.Q_exp  = breath.V_T/breath.t_exp/60;          % Expiratory flow rate  [cm3/s]

end
