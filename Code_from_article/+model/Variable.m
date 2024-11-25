classdef Variable < handle

    properties
        name
    end
    
    methods
        
        
        function obj = Variable(name)
            obj.name = name;
        end
        
    end
    
    methods(Abstract)
        [n] = NumVar(obj)
        [y0] = InitialState(obj)
        [m, vm] = ComputeMass(obj, t, y)
    end
    
end

