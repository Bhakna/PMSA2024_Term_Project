function [Kp_ri_out, Kp_po_out, Kp_adi_out, Kp_gut_out, Kp_sp_out, Kp_hep_out, Kp_lung_out, ver_vss] = GetUpdatedKpValues_human(par, V, Vss_obs)
% Calculates rat Kp-values based on initial Kp-values from simCYP and 
% observed Vss (Vss_obs)

% Human Kp-values predicted from SimCYP (method 2, Rodgers et al)
Kp_heart  = 3.732517; 
Kp_kidney = 3.307841; 
Kp_brain  = 1.715078; 
Kp_bone   = 1.305684; 
Kp_muscle = 3.494328; 
Kp_skin   = 1.796166; 

Kp_ri     = mean([Kp_heart, Kp_kidney, Kp_brain]);
Kp_po     = mean([Kp_bone, Kp_muscle, Kp_skin]);
Kp_adi    = 0.5275882; 
Kp_gut    = 3.511429;
Kp_sp     = 3.685111; 
Kp_hep    = 5.625239; 
Kp_lung   = par.kp_u*par.fu;
% Kp_nose   = par.Vu_lung*par.fu;

% Physiological values from Boger et al. -ref found in paper

% Calculation of Kp-factor

% a) Calculate the sum of V_i*Kp_i excluding the lung and the nose
sum_V_Kp   = sum([(V.adi*Kp_adi),  (V.g*Kp_gut), (V.h*Kp_hep), (V.ri*Kp_ri), (V.s*Kp_sp), (V.po*Kp_po)]);

% b) Sum up V_art and  V_vein
V_pl      = V.art + V.vein;

% c) Calculate the Kp-factor
kp_factor =  (Vss_obs - (V_pl + Kp_lung*V.lung))/sum_V_Kp;

% Update Kp-values
Kp_ri_out     = Kp_ri*kp_factor;
Kp_po_out     = Kp_po*kp_factor;
Kp_adi_out    = Kp_adi*kp_factor;
Kp_gut_out    = Kp_gut*kp_factor;
Kp_sp_out     = Kp_sp*kp_factor;
Kp_hep_out    = Kp_hep*kp_factor;
Kp_lung_out   = Kp_lung; 
% Kp_nose_out   = Kp_nose;

ver_vss = sum([(V.adi*Kp_adi_out),  (V.g*Kp_gut_out), (V.h*Kp_hep_out), (V.ri*Kp_ri_out), (V.s*Kp_sp_out), (V.po*Kp_po_out), ...
               (Kp_lung_out*V.lung), V_pl]);
           
end