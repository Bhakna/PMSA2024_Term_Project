classdef Dissolution

    properties
        var_part
        var_fluid
        var_gut
        
        transport_velocity
        D

        drug_D 
        dens  
        Cs 
        dr 
        fu_f 

        
    end
    
    methods
        function [obj] =  Dissolution( var_part, var_fluid, D, drug_data, init_transport_velocity, var_gut )
        
        obj.var_part  = var_part;
        obj.var_fluid = var_fluid;          
        obj.D         = D;
        
        obj.drug_D  = drug_data.D;
        obj.dens    = drug_data.dens;
        obj.Cs      = drug_data.Cs;
        obj.dr      = drug_data.dr;
        obj.fu_f    = drug_data.fu_f;

        if nargin == 6
        	obj.transport_velocity = init_transport_velocity;
            obj.var_gut = var_gut;        
        else
        	obj.transport_velocity = [];
            obj.var_gut = [];
        end
        
        end
        
    end
    
end

