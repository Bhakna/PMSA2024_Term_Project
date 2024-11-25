classdef Deposition
    % This class simulates the deposition at creation and returns particle
    % distributions via Normal
    %   The model is not guaranteed to give accurate results for very small
    %   radii
    
    properties
        dep_fraction_grid_a
        dep_fraction
    end
    
    methods
        % Input da_min, da_max should be in [dm]
        
        function obj = Deposition(data, breath, da_min, da_max)
            
            % Run deposition model over relevant particle sizes [um]
            
            obj.dep_fraction_grid_a  = logspace(log10(da_min*1e5),log10(da_max*1e5),1000);
            
            [obj.dep_fraction, ~] = deposition.DepositionPattern(data, breath, obj.dep_fraction_grid_a);
            
            % Convert grid back to [dm]
            obj.dep_fraction_grid_a  = obj.dep_fraction_grid_a * 1e-5; 
            
        end
        
        % Normally distributed particles
        function [r, p_0_u, p_0_l, ET_dose_nmol] = Normal(...
                        obj, mu, sigma, data, drug_data, ...
                            lung_dose_ug, x_u, x_l)


            % Setup grid with respect to da limits.
            unit_dens_g = 1000; % [g/dm3]
            shape_factor = 1;    % shape factor for spheres = 1
            da_to_d = 1/sqrt( (drug_data.dens_g/(unit_dens_g*shape_factor)) );
            d_min = da_to_d*obj.dep_fraction_grid_a(1) + eps;
            d_max = da_to_d*obj.dep_fraction_grid_a(end) - eps;
            
            r = logspace(log10(d_min/2),log10(d_max/2),1000);    
            
            if (r(1) > mu - 3*sigma) || (r(end) < mu + 3*sigma)
                error('mu +-3sigma should be in [da_min, da_max]')
            end
            
            % Particles are normally distributed with regards to diameter  
            y = normpdf(r, mu, sigma);
            d_a = 2*r/da_to_d;
            

            % Get deposition fractions on particle dist grid
            dep_frac = interp1(obj.dep_fraction_grid_a, obj.dep_fraction', d_a, 'spline')';

            % Calculate the amount of particles per radius deposited in each generation
            % Integration of y_m_dep wrt r (and not d_int!) is subsequently needed to 
            % get the amount deposited
            
            y_dep = dep_frac.*y;
            
            % Scale to lung deposited dose
            norm_factor = sum(trapz(r, y_dep, 2));
            y_m_dep = y_dep / norm_factor * lung_dose_ug*1e-6;
            

            % Next we convert the "generation" based formulation to a
            % continuous one

            % Add zero-radius
            y_m_dep = [interp1(r, y_m_dep', 0, 'PCHIP', 'extrap')', y_m_dep];
            r = [0, r];

            [N,M] = size(y_m_dep);

            % Compute mass
            mp = r.^3*4*pi/3 * drug_data.dens_g;

            MP = ones(N,1) * mp;

            y_part = y_m_dep./MP;

            % Set the mass corresponding to particles of size 0 to 0
            y_part(MP(:,:)==0) = 0;


            % Calculate amount of compound deposited [g]
            y_0      = y_part(2:end, :);        % Remove the ET-region
            y_0_L    = y_0./(data.L.value*ones(1,M));              % Number of particles deposited per radius [dm] per length [dm]


            p_0_u = zeros(length(x_u), M);
            p_0_l = zeros(length(x_l), M);


            for i = 1:M
                p_0_u(:,i) = interp1(data.L.grid, y_0_L(:,i), x_u, 'PCHIP');
                p_0_l(:,i) = interp1(data.L.grid, y_0_L(:,i), x_l, 'PCHIP');
            end

            % Normalize to retrieve the same dose as in the discrete model

            MP_u = (ones(size(p_0_u,1),1) * mp);
            MP_l = (ones(size(p_0_l,1),1) * mp);

            f = 1/( trapz(x_l, trapz(r, p_0_l.*MP_l, 2))*1e6 + trapz(x_u, trapz(r, p_0_u.*MP_u, 2))*1e6) * lung_dose_ug;

            p_0_l = p_0_l *f ;
            p_0_u = p_0_u *f ;

            % Compute extrathoracic (ET) and total lung deposition fraction
            lung_df = sum(trapz(r(2:end),y_dep(2:end,:),2)); % Total lung deposition fraction
            ET_df   = trapz(r(2:end),y_dep(1,:));
            ID_dose_ug   = lung_dose_ug/lung_df;
            ET_dose_ug   = ID_dose_ug*ET_df;
            ET_dose_nmol = 1e9*(ET_dose_ug*1e-6)/drug_data.MW;
            
        end
    end
end

