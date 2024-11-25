function [df, imax] = DepositionPattern(data, breath, da)
% Function for calculating human drug deposition accounting for bolus
% scaling
% Input: breath, a struct defined in BreathingParameters.m, and da [um]
% Output: deposition fractions (df), imax 

% A typical path model for human drug deposition accounting for bolus scaling

% Script for calculating deposition in the human lung after tidal breathing
% with possibility to assign: 1) different times for the inspiratory and the
% expiratory phases (t_in and t_exp, resp), and 2) a breathing pause.


%% Breathing conditions and anatomy

if (breath.V_D + breath.V_B) > breath.V_T
    error('The sum of the delay and bolus volume cannot exceed the tidal volume')
end

% Humans don't have any ET-compartments added to the anatomy
% data-set --> add an additional row at the beginning which
% subsequently is replaced by mouth deposition probabilities

% Add zero to represent oral impaction
N     = [0; data.N.value]; 
L     = [0; data.L.value]*10; % [cm]
D     = [0; data.D.value]*10; % [cm]

phi   = [0; data.phi.value];   % Angle 1, with gravity, [rad]
theta = [0; data.theta.value]; % Angle 2, branching angle, [rad]

V = L.*N*pi.*(D/2).^2;         % [cm3]

% Scale alveolar region
V_add  = data.V_lung_tot - sum(V);
V      = V + [0; data.alv_frac.value]*V_add;

%% Function for transforming volumes, diameters and length

[frac_bolus, frac_pause, L_scaled, D_scaled, imax, ~, V_2] = deposition.BolusScaling(breath, L, D, V);

% Choose input values
L = L_scaled;
D = D_scaled;

%% Particle input parameters

% Unit particle density
po = 1000; % 1 g/cm3 = 1000 kg/m3

%% Calculation of airflow (Q_i), velocity (V) and mean residence time (t_i)
    
[V_in,t_i_in,~]    = deposition.ComputeAirflow(D,N,breath.Q_in,V_2);     % Inspiratory
[V_exp,t_i_exp,~] = deposition.ComputeAirflow(D,N,breath.Q_exp,V_2);     % Expiratory

% Transformation of V to SI-units
V_in  = V_in./100; % from cm/s to m/s

V_exp = V_exp./100; % from cm/s to m/s

% Get deposition parameters 
[Cd,vg,ne,pa,Dmol] = deposition.ParametersA(da,po); 

%% Get deposition equation parameters for human (upright position) and calc deposition
% 1a) Inspiratory
[stk_in,re_in,eps_in, ~]      = deposition.ParametersB(da,po,ne,V_in,D,t_i_in,pa,vg,phi,Cd, breath.t_b);

% 1b) Expiratory
[stk_exp,re_exp,eps_exp, eps_b]   = deposition.ParametersB(da,po,ne,V_exp,D,t_i_exp,pa,vg,phi,Cd, breath.t_b);

% 2a) Inspiratory
[IMP_i_in, SED_in, DIF_i_in]    = deposition.DepositionProbability(stk_in,re_in,theta,eps_in, D, L,V_in,Dmol);

% 2b) Expiratory
[IMP_i_exp, SED_exp, DIF_i_exp] = deposition.DepositionProbability(stk_exp,re_exp,theta,eps_exp, D, L,V_exp,Dmol);

% 2c) Breath-hold
[IMP_i_b, SED_b, DIF_i_b]       = deposition.DepositionProbabilityBreathHold(eps_b, D, da, breath.t_b, phi);
    
% 3) Calculate oral deposition in man
[oral_in]  = deposition.ComputeOralDeposition(breath.Q_in, da, Dmol); % Cheng et al. 2003
[oral_out] = deposition.ComputeOralDeposition(breath.Q_exp, da, Dmol); % Cheng et al. 2003

% oral doesn't distinguish between different dep mechanism --> everything is placed on impaction 

% No oral deposition assumed during during breath-hold
IMP_i_in(1,:) = oral_in;
IMP_i_exp(1,:)= oral_out;

% Run deposition model
[ df, DEP_in, DEP_ex, DEP_b, df_in ] = deposition.ApplyFractions( IMP_i_in, SED_in, DIF_i_in,IMP_i_exp, SED_exp, DIF_i_exp, ...
                                                IMP_i_b, SED_b, DIF_i_b, V_2, imax, frac_bolus, frac_pause );                                            

% Pad df with zeros to cover unexposed generations
n     = length(N)-1;
[k,~] = size(df);                      
df(end+1:end+n-k+1,:) = 0;
