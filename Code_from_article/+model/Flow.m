classdef Flow
    %FLOW Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        var
        flow_rate
    end
    
    methods
        function [obj] = Flow(var0, var1, flow_rate)
            
            if ~isa(flow_rate, 'double')
                error('flow_rate must be a double!')
            end
            
        	obj.var{1,1} = VarName(var0);
            obj.var{2,1} = VarName(var1);
            obj.flow_rate = flow_rate;          
        end
    end
    
end

function [name] = VarName(var)
    switch class(var)
        case 'char'
            name = var;
        case 'model.Concentration'
            name = var.name;
        otherwise
            error('Attempted to add unsupported variable to flow!')
    end

end