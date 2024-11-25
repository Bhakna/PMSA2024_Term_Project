function [ df, DEP_in, DEP_ex, DEP_bb, df_in ] = ApplyFractions( IMP_in, SED_in, DIF_in,IMP_exp, SED_exp, DIF_exp,IMP_b, SED_b, DIF_b, V, imax, frac_bolus, frac_pause )
% Calculation of deposition fractions (df) when bolus-dosing is applied 
% Different df for inspiration and expiration

% Calculation of resulting probability P
% a) Inspiratory
P_in = 1-(1-IMP_in).*(1-SED_in).*(1-DIF_in);

% b) Expiratory
P_exp = 1-(1-IMP_exp).*(1-SED_exp).*(1-DIF_exp);

% c) Breath-hold
P_b = 1-(1-IMP_b).*(1-SED_b).*(1-DIF_b);

N_g = imax;         % Nr of generations
N_p = size(P_in,2); % Nr of particle sizes

% Augment P to compute f, set the first row to 0
% --> first row in fi should equal 1 (nothing has deposited when entering
% the nose)

P_hat = [zeros(1,N_p); P_in(1:imax,:)]; 

f = zeros(N_g+1,N_p);

for i = 1:(N_g+1)
    f(i,:) = prod([ones(1,N_p);1-P_hat(1:i,:)]); 
end

% Clone V into a matrix
V = V(1:imax)*ones(1,N_p);
Frac_pause = frac_pause(1:imax)*ones(1,N_p);

% Exhalation
x = zeros(N_g, N_p);

% Mass left at bottom
y = zeros(N_g, N_p);

% Create beta for breath-hold
beta = f(3:end,:).*Frac_pause(2:end,:);

% x defines all that flows, the probability for deposition is taken into
% account afterwards in DEP_ex: DEP_ex(i)=P(i)*x(i)

for i = N_g-1:-1:1 

    x(i,:)   = (1-P_exp(i+1,:)).*x(i+1,:) + (1-P_b(i+1,:)).*beta(i,:); %*(1-P_b(i+1,:))
    y(i+1,:) = P_b(i+1,:).*beta(i,:);
    
end

DEP_ex = x.*P_exp(1:imax,:);

% Inhalation
% Create vector of cumulative volumes (take away V that stays in each
% generation)

sum_V=zeros(length(N_g),1);
for i = 1:N_g
    V2=flipud(V(:,1));
    sum_V(i)=sum(V2(1:(N_g+1-i)));
end

DEP_in = zeros(N_g, N_p);
for j = 1:N_p

    DEP_in(:,j) = f(1:end-1,j).*P_in(1:imax,j).*frac_bolus(1:N_g);

end

% Breath-hold
DEP_bb = zeros(N_g, N_p);
for j = 1:N_p % for each particle size
    DEP_bb(:,j) = f(1:end-1,j).*(1-P_in(1:imax,j)).*P_b(1:imax,j).*V(:,1); % add the V as well
end

% Result including breath-hold:
df = DEP_in + DEP_ex + DEP_bb;
df_in = DEP_in;

