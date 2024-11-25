classdef Constant

    properties
        value
        grid
    end
    
    methods
        function obj = Constant(var, grid)
            switch class(var)
                case 'double'
                    obj.value = var;
                    if nargin < 2
                        obj.grid = [];
                    else
                        obj.grid  = grid;
                    end
                case 'model.Constant'
                    obj = var;
                otherwise
                    error('Unsupported input in constant')
            end
        end

        function [new_obj] = ToGrid(obj, grid)
            if ( ~isempty(grid) && ~isempty(obj.grid) )
                new_obj = obj;
                new_obj.value = interp1(obj.grid, obj.value, grid(:),'PCHIP','extrap');
                new_obj.grid = grid(:);
            end
        end

        function [value] = Integral(obj, grid)
            value = trapz(grid, obj.OnGrid(grid));
        end

        function [value] = AsVolumes(obj,grid)
           dx = [diff(grid); 0]/2 + [0;diff(grid)]/2;
           value = obj.OnGrid(grid).*dx;
        end
        
        function [value] = OnGrid(obj, grid)
            if ( ~isempty(grid) && ~isempty(obj.grid) )
                value = interp1(obj.grid, obj.value, grid,'PCHIP','extrap');
            elseif ( isempty(grid) && isempty(grid) )
                value = obj.value;
            else
                error('!!!')
            end
        end

        function [val1, val2, grid] = SameGrid(obj1, obj2 )
           
            [val1, grid1] = GetValueAndGrid(obj1);
            [val2, grid2] = GetValueAndGrid(obj2);
           
            if isa( obj1, 'double') 
                grid = grid2;
            elseif isa( obj2, 'double')
                grid = grid1;
            else
                if isequal(grid1, grid2)
                    grid = grid1;
                else
                    error('Constants on different grids not supported!')
                end
            end
        end

        
        function r = uplus(obj1)
             r = obj1;
        end
        function r = uminus(obj1)
            [val, g] = GetValueAndGrid(obj1);
             r = model.Constant( -val, g );
        end
        function r = plus(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1+val2, g );
        end
        function r = minus(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1-val2, g );
        end
        function r = mtimes(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1*val2, g );
        end
        function r = mrdivide(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1/val2, g );
        end
        function r = mpower(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1^val2, g );
        end
        function r = times(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1.*val2, g );
        end
        function r = rdivide(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1./val2, g );
        end
        function r = power(obj1,obj2)
             [val1, val2, g] = SameGrid( obj1, obj2 );
             r = model.Constant( val1.^val2, g );
        end
    end
    
end

function [val, grid] = GetValueAndGrid(object)
    if ( isa( object, 'double') )
        val = object;
        grid = [];
    elseif ( isa( object, 'model.Constant') )
        val = object.value;
        grid = object.grid;
    else
        error('Unsupported type!')
    end
end