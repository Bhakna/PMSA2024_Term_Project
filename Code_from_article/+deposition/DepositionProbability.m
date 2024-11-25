function [DEP_i, SED_i, DIF_i]=DepositionProbability(stk,re,theta,eps, D_i, L_i,v_i,Dmol)
% Calculation of deposition probability for the airway tree during inhalation 
% pr exhalation

D_i     = D_i.*10^-2; % Transformation from cm to m
L_i     = L_i*10^-2;  % Transformation from cm to m

[n, m]=size(stk); % n = number of gen; m = number of particle sizes


DEP_i=zeros(n,m);
SED_i=zeros(n,m);
DIF_i=zeros(n,m);
delta=zeros(n,m);


% Impaction from Yu and Diu 1982
for j=1:m
    for i=1:n
        
        theta_temp = L_i(i)/(4*D_i(i));
        DEP_i(i,j) = 0.768*theta_temp*stk(i,j);                       

    end
end

% Sedimentation from Thomas 1959
for j=1:m
    for i=2:n

        SED_i(i,j) = 2/pi*(2*eps(i,j)*(1-eps(i,j)^(2/3))^(1/2)...
                     -(eps(i,j)^(1/3))*(1-eps(i,j)^(2/3))^(1/2)...
                     +asin(eps(i,j)^(1/3)));

    end
end

% Remove any complex numbers
SED_i = real(SED_i);

% Diffusion from Ingham 1975
% Differences for turbulent and laminar flow according to Yu and Diu 1982

for j=1:m
    for i=2:n
        delta(i,j)=Dmol(j)*L_i(i)/(v_i(i)*(D_i(i)^2));
        
        if re(i,j) > 2000 % turbulent flow
            DIF_i(i,j)=4*(delta(i,j)^0.5) * (1 - 0.444*(delta(i,j)^0.5)); % the eq. continues with a (...) in Yu and Diu 1982
        else              % laminar flow
            DIF_i(i,j)=1-0.819*exp(-14.63*delta(i,j))-0.0976*exp(-89.22*delta(i,j))-0.0325*exp(-228*delta(i,j))...
                -0.0509*exp(-125.9*delta(i,j)^(2/3));
        end
    end
end
