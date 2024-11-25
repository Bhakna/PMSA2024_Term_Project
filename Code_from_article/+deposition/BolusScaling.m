function [frac_bolus, frac_pause, L_scaled, D_scaled, imax, V_cum_scaled, V_scaled] = BolusScaling(breath, L, D, V)
% Function that returns scaled airway geometries and volume fractions according 
% to the bolus principle described in:
% Bondesson et al., J Aerosol Med Pulm,  vol 18(1), pp. 23-33, 2005
% Input: breath, a struct defined in BreathingParameters.m

TLC     = sum(V);  % Total lung capacity, TLC [mL]

% 1) Scale volume to desired FRC
scale_init = breath.FRC/TLC;
Vi_init    = V*scale_init;

% 2) Check for unrealistic values
if (breath.V_D + breath. V_B) > breath.V_T
    error('The sum of the delay and bolus volume cannot exceed the tidal volume')
end

% Create a vector with cumulative volumes
cum_V = cumsum(Vi_init);
cum_Vi_init = sum(Vi_init)*ones(length(Vi_init),1)- [0;cum_V(1:end-1)];

% Scaling of the front (scale_f)
scale_f = ones(length(Vi_init),1);

for i = 1:length(Vi_init)
    
    scale_f(i) = (sum(Vi_init) + breath.V_D)/cum_Vi_init(i);
    
end

% Scaling of the tail (scale_t)
scale_t = ones(length(Vi_init),1);

for i = 1:length(Vi_init) 
    
    scale_t(i) = (sum(Vi_init) + breath.V_D + breath.V_B)/sum(Vi_init(i+1:end));
    
    % Special treatment for the last generation
    if i == length(Vi_init) 
        scale_t(i) = (sum(Vi_init) + breath.V_D + breath.V_B)/sum(Vi_init(i:end));
    end
    
end

% Check if the V_T is insufficient to wash out the bolus tail
checking = zeros(length(Vi_init)-1,1);

for i = 1:length(Vi_init)-1
    
    checking(i) = breath.V_T < ((cum_V(i) + breath.V_D + breath.V_B)/(1 - cum_V(i)/sum(Vi_init))); % This formula accounts for scaling

end


% Calculate the first generation where the bolus is not washed out       
[val_max, i_wash] = max(checking);      

if val_max == 1 && i_wash == 1
    i_wash = length(Vi_init);
end

% Different calculation of scale_t between i_wash and imax
scale_t(i_wash:end) = (sum(Vi_init) + breath.V_T)/sum(Vi_init);

% Calculate average scale factor (f_ave) over front and tail
f_ave = (scale_f + scale_t)/2;

V_scaled = f_ave(1:end).*Vi_init(1:end);
V_cum_scaled = cumsum(V_scaled);

% Calculate the last generation ventilated by the bolus (imax)
[vaLmax, imax] = max(breath.V_T - breath.V_D < (V_cum_scaled(1:end-1)));

if vaLmax == 0 && imax == 1
    imax = length(Vi_init);
end

% 3) Adapt volume of imax such that sum_Vi=V_T
V_scaled(imax) = (breath.V_T - breath.V_D) -V_cum_scaled(imax-1);

% 4) Set any generations after imax to the corresponding initial value
if imax < length(Vi_init)
    V_scaled(imax+1:end) = Vi_init(imax+1:end);
end

% 5) Update V_cum_scaled and airway dimensions
V_cum_scaled = cumsum(V_scaled);

L_scaled   = L.*f_ave.^(1/3);
D_scaled   = D.*f_ave.^(1/3);


% 6a) Fraction of bolus reaching generation i, bolus_frac
frac_bolus             = zeros(length(V),1);
frac_bolus(1:i_wash-1) = ones(i_wash-1,1);

VF_cum = zeros(length(V),1);
for i = i_wash:imax
    VF_cum(i) = (breath.V_T - breath.V_D - V_cum_scaled(i-1))/breath.V_B;
end
VF_cum(VF_cum>1) = 1;

% 6b) Special treatment for generations where the bolus does not completely
% pass
frac_bolus(i_wash:imax) = VF_cum(i_wash:imax);

% 7) Volume fraction at pause (frac_pause)
frac_pause = zeros(length(frac_bolus),1);
for i = 2:imax
    
    if i < imax
        frac_pause(i) = frac_bolus(i) - frac_bolus(i+1);
    elseif i == imax
        frac_pause(i) = 1 - sum(frac_pause);
    end
    
end


end

