function [ data ] = Scale( data, factor )

fields = fieldnames(data);

for i=1:numel(fields)
    data.(fields{i}) = factor * data.(fields{i});
end

end

