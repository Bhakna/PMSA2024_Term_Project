function [stk,re,eps, eps_b]=ParametersB(da,po,ne,v_i,D_i,t_i,pa,vg,phi,Cd, t_b)
% Calculation of deposition parameters: eps, reynolds number, stokes
% number and eps

% eps is calculated separately for breath-hold depending on breath-hold time (t_b)

% v_i already written in the correct units
da  = da.*10^-6;  % Transformation from um to m
D_i = D_i.*10^-2; % Transformation from cm to m

m = length(da);
n = length(v_i);

stk   = zeros(n,m); % Stokes number
re    = zeros(n,m); % Reynolds number
eps   = zeros(n,m); % eps during inhalation and exhalation
eps_b = zeros(n,m); % eps_b during breath-hold

% According to Yu and Diu  1982, all phi are set to pi/4

phi = (pi/4).*ones(size(phi)); % sin(0.9033).*ones(size(phi));  

for j=1:m
    for i=1:n

            stk(i,j) = po*(da(j)^2)*v_i(i)*Cd(j)/(9*ne*D_i(i)); % With cunningham, Zhang et al. 1997
            re(i,j)  = pa*D_i(i)*v_i(i)/ne;

            eps(i,j)   = 3*vg(j)*t_i(i)*cos(phi(i))/(4*D_i(i)); % upright position
            eps_b(i,j) = 3*vg(j)*t_b*cos(phi(i))/(4*D_i(i));    % upright position

    end
end

end

% eps = par in sedimentation equation
% re  = Reynolds number
% stk = stokes number