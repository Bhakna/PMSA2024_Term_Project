function [ flows ] = DualFlow( var0, var1, q0, q1 )

flows = [ model.Flow(var0, var1,  q0);
          model.Flow(var1, var0,  q1)];

end

