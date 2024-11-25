function [ data ] = loadData( files )

if isa( files, 'char' )
    files = {files};
elseif  ~isa( files, 'cell' )
    error('Unsupported file specification.')
end

data = [];


for i = 1:length(files)

    temp = importdata(['+phys/' files{i}]);

    if size(temp.textdata,1) == 1
        
        if size(temp.data,1) == 1
            for j = 1:length(temp.colheaders)
                data.(temp.colheaders{j}) = temp.data(j);
            end 
        else
            grid = temp.data(:,1);
            for j = 2:length(temp.colheaders)
                data.(temp.colheaders{j}) = model.Constant( temp.data(:,j), grid );
            end 
        end
    else
        for k = 2:size(temp.textdata, 2)
            for j = 2:size(temp.textdata, 1)
            
                data.(temp.textdata{1,k}).(temp.textdata{j,1}) = temp.data(j-1, k-1);
            end
        end
    end
    
end


end

