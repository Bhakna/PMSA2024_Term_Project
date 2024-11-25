classdef Concentration < model.Variable

    properties
        vol
        y0
        grid
        dx
    end
    
    methods
        function obj = Concentration(name, vol, grid, y0)
            
            obj@model.Variable(name);
            
            switch class( vol )
                case 'model.Constant'
                    error('Make sure to imput raw double data and not model.Constant!')
                case 'double'
                    obj.vol = vol;
                otherwise
                    error('Unsupported volume input!')
            end

            
            switch nargin
                case 3
                    obj.y0 = zeros(size(grid(:)));
                    obj.grid = grid(:);
                case 4
                    obj.y0 = y0;
                    obj.grid = grid(:);
                otherwise
                    obj.y0 = 0;
                    obj.grid = [];
            end
            
            
            if ~isempty(vol) && ~isequal(length(vol),length(vol))
                error('Volume and grid no same length!')
            end

            if length(obj.grid) < 1
                obj.dx = 1;
            else
                obj.dx = [0;diff(obj.grid)/2] + [diff(obj.grid)/2;0];
            end
        end
        
        function [dx] = Length(obj)
            dx = obj.dx;
        end
        
        function [volume] = Volume(obj)
            volume = sum(obj.vol);
        end

        function [] = ScaleToVolume(obj, volume)
            obj.vol = obj.vol / sum(obj.vol).*volume;
        end
        
        function [n] = NumVar(obj)
            n = length(obj.vol);
        end
        function [mass, vmass] = ComputeMass(obj, t, y)
            vmass = obj.vol'.*y;
            mass  = sum(vmass,2);
        end
                
        function [y0] = InitialState(obj)
            y0 = obj.y0;
        end
    end
    
end

