function [v_i,t_i,Q_i] = ComputeAirflow(D_i_scaled,N_i,Q,V_i)
% Calculation of air velocity (v_i) and residence time (t_i) in each tube
% Prior to this, calculate the flow Q_i for each generation i

% 1) Calculation of the flow in each airway generation i (Q_i)
Q_i = Q./N_i; % [cm3/s]

% 2) Calculation of the cross-sectional area in each airway generation i (A_i)
A_i = pi.*(D_i_scaled./2).^2; % [cm2] A=r^2*pi

% 3) Calculation of the velocity in each airway generation (v_i)
v_i = Q_i./A_i; % [cm/s]

% 4) Calculation of the residence time in each airway generation (t_i)
t_i    = V_i./N_i./Q_i; % [s]
t_i(1) = 0;

% Due to additional alveolar volume, t_i is historically calculated in this way
% rather than t_i = L_i/v_i
% As v_i = Q_i./A_i --> t_i = L_i/v_i = L_i*A_i/Q_i = V_i/N_i/Q_i
end

