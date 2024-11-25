classdef Particles < model.Variable

    properties
        x_grid
        z_grid
        drug_dens
        p0
        m
        dissolution
        
    end
    
    methods
        function [obj] = Particles(name, x_grid, z_grid, drug_dens, p0 )
            obj@model.Variable(name);
            
            obj.x_grid = x_grid;
            obj.z_grid = z_grid;
            obj.drug_dens = drug_dens;        
            obj.p0 = p0;
            
            obj.dissolution = [];
            
            % Precompute mass
            
            m = zeros(size(p0));
            for j = 1:length(z_grid)
                z = max(0, z_grid - z_grid(j));
                m(:,j) = trapz(z, p0.*z*2/3*pi, 2);
            end
            
            obj.m = griddedInterpolant({x_grid, z_grid}, m, 'spline', 'nearest');


        end
        
        function [] = SetDissolution(obj, item)
        	obj.dissolution = item;
        end
        function [n] = NumVar(obj)
            n = length(obj.x_grid);
        end
        
        function [mass, vmass] = ComputeMass(obj, t, y)
           
            if isempty(obj.dissolution)
            	error('Model was not compiled yet!')
            else
            
                N = size(y,2);

                if obj.dissolution.pre.simple
                    vmass = obj.m(ones(length(t),1)*obj.x_grid',-y)*obj.drug_dens;
                else
                    source = obj.dissolution.pre.x0_interp(t*ones(1,N), ones(length(t),1)*(1:N));
                    nu = obj.dissolution.pre.nu_interp(t*ones(1,N), ones(length(t),1)*(1:N));
                    vmass = obj.m(source,-y)*obj.drug_dens.*nu;
                end
                mass = trapz(obj.dissolution.var_p.x_grid, vmass, 2);


            end
        end
        
        function [y0] = InitialState(obj)
            y0 = zeros(obj.NumVar(), 1);
        end
    end
    
end

