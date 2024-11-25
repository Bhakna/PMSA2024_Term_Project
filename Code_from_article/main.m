close all
clear variables


case_study = @ cases.CaseStudyOne; % @ cases.CaseStudyTwo; % @ cases.CaseStudyThree;

anatomy = 'lung_weibel.csv'; % 'lung_yeh_schum.csv';


% Load physiology
files = {anatomy, 'pbpk.csv', 'scalars.csv'};
[ data ] = phys.loadData( files );
data.V   = phys.Scale( data.V_f, data.BW );
data.Q   = phys.Scale( data.f, data.Q_co );

% Load drug data
[ drug_data ] = drug.loadData(data.V_f, data.BW, data.al_ind, data.gen.value(end) );

% Setup lung grid
N_upper = 100;
N_lower = 100;

x_b = [0, data.L.grid(data.al_ind), data.L.grid(end)];

x_upper = linspace(x_b(1), x_b(2), N_upper)';
x_lower = linspace(x_b(2), x_b(3), N_lower)';

% Setup deposition model
da_min = 0.00015*1e-5; % Minimum diameter [dm]
da_max = 20*1e-5;     % Maximum diameter [dm]

breath = deposition.BreathingParameters();
dep = deposition.Deposition(data, breath, da_min, da_max);


% Run case study
for i = 1:3

    % Load case study data
    [mu, sigma, drug_data.Cs, T, init_transport_velocity,...
                lung_dose_ug, plotFunction] = case_study(i);


    % Create particle distribution
    [r, p_0_u, p_0_l, ET_dose_nmol] = dep.Normal(mu, sigma, data,...
                             drug_data, lung_dose_ug, x_upper, x_lower);


    % Create model
    m = model.Model();
    CreateModel

    % Compile and run
    c = m.Compile(T);
    [ t, y ] = simulate.RunSimulation(c, T);

    % Post-simulation (down sample t/y for speed if needed)
    [mass, vmass] = m.ComputeMass( t, y );
    [conc] = m.GetConcentrations( t, y );

    plotFunction(conc, mass, vmass, data, drug_data, ...
                    lung_dose_ug, x_lower, x_upper, t, i);


end
