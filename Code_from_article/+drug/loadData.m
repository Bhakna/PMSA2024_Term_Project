function [data] = loadData(volumes, BW, al_ind, N_gen)
% Function defining drug-specific input drugameters for salbutamol in humans
% PK-drugameters

data.CL      = 1;                 % Plasma CL [L/h/kg]
data.fu      = 0.75;              % Fraction unbound in plasma 
data.b_p     = 1;                 % Blood-plasma ratio
data.Cl      = data.CL/data.b_p;  % Blood-clearance
data.Cl      = data.Cl*BW;        % Transformation to L/h
data.fu_f    = 1;                 % Free fraction in the ELF
data.F       = 0.05;              % Oral bioavailability
data.ka      = 0.6;               % Oral absorption rate constant [h-1]

% Phys-chem
data.Papp    = 1.5;               % Apparent permeability [10-6 cm/s]
data.P_app   = data.Papp*10^(-6);  % unit [10-6 cm/s]
data.P       = 3600/10*data.P_app; % unit transform from cm/s to dm/h
data.MW      = 250;                % Molecular weight [g/mol]

[P_eff]      = drug.CalcPeffFromPapp(data.MW, data.P_app); % Effective permeability Sjögren et al. 2013
data.P_eff   = P_eff;
data.P_eff_unit = 3600/10*P_eff;       % unit transform from cm/s to dm/h
data.P       = data.P_eff_unit;        % Update the permeability

% Diffusion coefficient 
visc_factor  = 1;                                    % Viscosity in ELF relative to saline
[D_ES]       = drug.StokesEinstein(data.MW, visc_factor); % Get diffusion coefficient using Stokes-Einstein eq.
data.D       = D_ES;

data.kp_u    = 6.5; % Vu,lung
Vss_obs      = 2;      % Steady-state volume of distribution [L/kg]
[data.Kp.ri, data.Kp.po, data.Kp.adi, data.Kp.gut, data.Kp.sp, data.Kp.hep, Kp_lung_out, ~] = drug.GetUpdatedKpValues_human(data, volumes, Vss_obs);

data.k_p = [data.Kp.sp, data.Kp.ri, data.Kp.po, data.Kp.adi, data.Kp.hep,data.Kp.gut, Kp_lung_out];      % {'sp','ri','po','adi','hep','gut', 'lung'};

% Kp,u-values 
data.kp_u_all = data.k_p./data.fu;

% Kp-values 
data.Kp.lung  = data.kp_u*data.fu;

% Kp,u-values 
data.Kp_u.ri   = data.Kp.ri/data.fu;
data.Kp_u.po   = data.Kp.po/data.fu;
data.Kp_u.adi  = data.Kp.adi/data.fu;
data.Kp_u.gut  = data.Kp.gut/data.fu;
data.Kp_u.sp   = data.Kp.sp/data.fu;
data.Kp_u.hep  = data.Kp.hep/data.fu;

% Formulation-specific properties
data.dens_g   = 1000;                    % particle density [g/dm3]
data.dens     = data.dens_g/data.MW*1e9; % particle density [nmol/dm3]

tb_ind        = 1:al_ind;
al_ind        = al_ind+1:N_gen; 

% Scaling of the permeability
scale_al       = 10;             % Factor used for scaling the permeability in the al-region 10 --> 10-fold higher P
scale_tb       = 1;              % No scaling applied for the tb-region

data.scaleP        = ones(N_gen,1); % No scaling applied for the tb-region
data.scaleP(tb_ind)= scale_tb*ones(length(tb_ind),1);
data.scaleP(al_ind)= scale_al*ones(length(al_ind),1);

% Update data.P
data.P             = 2*data.P*data.scaleP;

% Other stuff
data.dr = 1;    % Particle layer as fraction of radius

end