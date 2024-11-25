function [ t, y ] = RunSimulation( comp, T )


options = odeset('Stats','on','JPattern', jpattern(comp));

tic
[t,y] = ode23t(@(t,y) odefunc(t,y,comp), [0,T], comp.y0, options);
toc

end

function [ dydt ] = odefunc(t, y, comp)

    dydt = comp.A * y; % Linear equations

    for i = 1:length(comp.diss)

        diss = comp.diss{i};

        [N_x, N_a] = size(diss.var_p.p0);

         % Begin with particle dissolution rate   
        dydt(diss.p_ind) = -2*diss.con.drug_D/diss.con.dens/diss.con.dr...
                *(diss.con.Cs - diss.con.fu_f *y(diss.f_ind));

        % Next is mass transfer
        if diss.pre.simple 
            
            y_int = zeros(N_x,1);
            for j = 1:N_x
               y_int(j,1) = diss.pre.p0interp{j}(-y(diss.p_ind(j)));
            end
            nu = 1;        

        else

            particle_source = diss.pre.x0_interp(t*ones(N_x,1), (1:N_x)');

            y_int = diss.pre.p0interp_2d(particle_source, -y(diss.p_ind));
            nu = diss.pre.nu_interp(t*ones(N_x,1), (1:N_x)');
        end

 
        dydt(diss.f_ind) = dydt(diss.f_ind) - ...
               diss.con.dens*2/3*pi*dydt(diss.p_ind).*y_int.*nu./diss.var_f.vol.*diss.var_f.dx;



        if ~diss.pre.simple

            % Current position of particles on p_grid        
            dadx = diff([y(diss.p_ind);y(diss.p_ind(end))])./diss.pre.dx;       

            dydt(diss.p_ind) = dydt(diss.p_ind) - dadx.*diss.pre.dxdt;


            % Now, mass transfer into gut        
            z = diss.var_p.z_grid' + y(diss.p_ind(1));

            [~, j] = max( z > 0 );

            if j > 1    % Unless particles are not completly disolved

                p_source_gut = particle_source(1);
                p0 = interp1(diss.var_p.x_grid, diss.var_p.p0, p_source_gut)';


                p0(j-1) = (0 - z(j-1))/(z(j) - z(j-1)) * ( p0(j-1) - p0(j) ) + p0(j-1);
                z(1:j-1) = 0;

                v0 = -diss.pre.dxdt(1);
                d0 = nu(1); 
                dydt(diss.pre.gut_ind) = dydt(diss.pre.gut_ind) + trapz(z, p0.*z * 2/3*pi * diss.con.dens)*d0*v0;
            end

        end

    end
    
end



function [pattern] = jpattern(comp)

    pattern = comp.A~=0;
    for i = 1:length(comp.diss)

        diss = comp.diss{i};

        for j = 1:length(diss.f_ind)
            pattern(diss.f_ind(j), diss.f_ind(j)) = 1;  
            pattern(diss.f_ind(j), diss.p_ind(j)) = 1;
            pattern(diss.p_ind(j), diss.f_ind(j)) = 1;
        end
        
        if ~diss.pre.simple
            for j = 1:length(diss.p_ind)
                pattern(diss.p_ind(j), diss.p_ind(j)) = 1;
                if j < length(diss.p_ind)
                    pattern(diss.p_ind(j), diss.p_ind(j+1)) = 1;
                end
            end
            pattern(diss.pre.gut_ind, diss.p_ind(1)) = 1;

        end

    end
end