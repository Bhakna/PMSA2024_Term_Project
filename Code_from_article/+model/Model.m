classdef Model < handle
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)

    end  
    properties (GetAccess = private)
        variables
        flows
        dissolutions
        
        comp
    end
    
    methods
        function obj = Model()
            obj.variables    = containers.Map();
            obj.flows        = model.Flow.empty;
            obj.dissolutions = model.Dissolution.empty;
            
            obj.comp = [];
        end
        
        function [] = Add(obj, item)

            if length(item) > 1
                for i = 1:length(item)
                    obj.Add(item(i));
                end
            else
                switch class(item)
                    case 'model.Concentration'
                        obj.variables(item.name) = item;
                    case 'model.Particles'
                        obj.variables(item.name) = item;
                    case 'model.Flow'
                        obj.flows(end+1) = item;
                    case 'model.Dissolution'
                        obj.dissolutions(end+1) = item;
                    otherwise
                        error('Attempted to add unsupported item to model!')
                end
            end
        end
        
        
        function [c] = Compile(obj, T)

            obj.CompileVariables();
            obj.CompileFlows();
            obj.CompileDissolutions(T);
            
            c = obj.comp;

        end
       
        function [mass, vmass] = ComputeMass( obj, t, y )
            labels = keys(obj.variables);
            
            mass = [];
            mass.total = zeros(size(t));
            
            vmass = [];
            
            for var = labels
               [mass.(var{1}),vmass.(var{1})] = obj.variables( var{1} )...
                            .ComputeMass( t, y(:,obj.comp.varmap(var{1})));
                        
                mass.total = mass.total + mass.(var{1});
            end
            
        end
        
        function [conc] = GetConcentrations( obj, t, y )
            labels = keys(obj.variables);
            
            conc = [];
            
            for var = labels
               conc.(var{1}) = y(:,obj.comp.varmap(var{1}));
            end
            
        end
        
        function [] = CompileVariables(obj)

            obj.comp.varmap = containers.Map();
            
            N = 0;
            obj.comp.y0 = [];
            
            for k = keys(obj.variables)
                var = obj.variables(k{1});
                n   = var.NumVar();
                
                obj.comp.varmap(var.name) = (N + 1: N + n)';
                N = N + n;
                
                obj.comp.y0(end+1:end+n) = var.InitialState();
            end
            
            obj.comp.NumVar = N;
            
            

        end
        function [] = CompileFlows(obj)
           
            obj.comp.A = zeros(obj.comp.NumVar);


                for j = 1:length(obj.flows)
                    
                    d = obj.ExtractFlowData(obj.flows(j));
                    
                    n1 = length(d.ind{1});
                    n2 = length(d.ind{2});
                    for i = 1 : max( n1, n2 )
                        source_i = min(i, length(d.ind{1}));
                        sink_i   = min(i, length(d.ind{2}));
                        
                        source = d.ind{1}(source_i);
                        sink   = d.ind{2}(sink_i); 
                        
                        obj.comp.A(source, source) = obj.comp.A(source, source)...
                                - d.rate(i)/d.vol{1}(source_i);
                            
                        obj.comp.A(sink,   source) = obj.comp.A(sink,   source)...
                                + d.rate(i)/d.vol{2}(sink_i);
                    end

                end
        end
       
        
        function [data] = ExtractFlowData(obj, flow)

            data = [];
            data.ind = {};
            data.vol = {};
            data.grid = {};
            data.dx = {};
                        
            N = zeros(3,1);
            for i = 1:2
                data.ind{i} = obj.comp.varmap(flow.var{i});
                data.grid{i} = obj.variables(flow.var{i}).grid;
                data.vol{i} = obj.variables(flow.var{i}).vol;
                N(i) = length(data.grid{i});
            end
            
            N(3) = length(flow.flow_rate);
            
            if (sum( (N == 0) | (N == max(N)) ) < 3)
                error(['Flow from ', flow.var{1} ' to ', flow.var{2},...
                        ' are on different grids!'])
            end

            data.rate = flow.flow_rate;
            
        end

        function [] = CompileDissolutions(obj, T)
            
            obj.comp.diss = [];
            
            for i = 1:numel(obj.dissolutions)
                diss = [];
                diss.con = obj.dissolutions(i);
                diss.var_p = obj.variables(diss.con.var_part);
                diss.var_f = obj.variables(diss.con.var_fluid);
                diss.p_ind = obj.comp.varmap(diss.con.var_part);
                diss.f_ind = obj.comp.varmap(diss.con.var_fluid);
                diss.con.D = diss.con.D.ToGrid(diss.var_p.x_grid);
                
                if ~isempty(diss.con.transport_velocity)
                    x_grid = diss.var_p.x_grid;
                    y0 = [x_grid; ones(size(x_grid))];
        
                    D = diss.con.D.OnGrid(x_grid);
                    cross_area = pi.*(D/2).^2;

                    % Transport velocity on x_grid
                    diss.pre.dxdt = - diss.con.transport_velocity*(cross_area - cross_area(end))/(cross_area(1) - cross_area(end));

                    % Numeric derivative
                    eps = 1e-4;
                    dDdx = @(x) (diss.con.D.OnGrid(x + eps) - diss.con.D.OnGrid(x - eps))/eps/2;

                    diss.pre.d2xdtdx = - diss.con.transport_velocity*pi*D.*dDdx(x_grid)/2/(cross_area(1) - cross_area(end));
                    % dissolution.d2xdtdx = 0*dissolution.d2xdtdx;
                    
                    options = odeset();
                    [t,y] = ode45(@(t,y) TransportOde(t,y,diss), [0,T], y0, options);

                    N = length(x_grid);
                    
                    diss.pre.t = t;
                    diss.pre.x0 = y(:,1:N);
                    diss.pre.nu = y(:,N+1:end);
                    diss.pre.x0_interp = griddedInterpolant({t,(1:N)},diss.pre.x0);
                    diss.pre.nu_interp = griddedInterpolant({t,(1:N)},diss.pre.nu);
                    diss.pre.dx = diff(diss.var_p.x_grid);
                    diss.pre.dx(end+1) = diss.pre.dx(end);
                    diss.pre.gut_ind = obj.comp.varmap(diss.con.var_gut);
                    diss.pre.simple = false;
                else
                    diss.pre = [];
                    diss.pre.simple = true;
                end
            
                
                x_I_0 = cumtrapz(diss.var_p.z_grid, diss.var_p.p0, 2);
                x_I_0 = x_I_0(:,end) - x_I_0;
                
                for j = 1:size(diss.var_p.p0,1)
                    diss.pre.p0interp{j} = griddedInterpolant(diss.var_p.z_grid, x_I_0(j,:), 'pchip', 'nearest');
                end
                diss.pre.p0interp_2d = griddedInterpolant({diss.var_p.x_grid, diss.var_p.z_grid}, x_I_0, 'spline', 'nearest');

                diss.var_p.SetDissolution(diss);
                obj.comp.diss{i} = diss;
            end
            

        end
        
    end
    
end


function [dydt] = TransportOde(t,y,diss)



    [N_x, N_a] = size(diss.var_p.p0);
    
    dydt = zeros(2*N_x,1);
      
    % Transport velocity at slices
    dydt(1:N_x) = -interp1(diss.var_p.x_grid, diss.pre.dxdt, y(1:N_x), 'PCHIP');

    % d2xdtdx is the derivative of transport velocity
    dydt(N_x+1:2*N_x) = -y(N_x+1:2*N_x) .* interp1(diss.var_p.x_grid, diss.pre.d2xdtdx, y(1:N_x), 'PCHIP');

end
