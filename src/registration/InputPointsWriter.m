function InputPointsWriter( filename,points,varargin )
%INPUTPOINTSWRITER Write a file with points that can be read by
%transformix.


% Read inputs
p = inputParser;
addRequired(p,'filename')
addRequired(p,'points')
addParameter(p,'type','index',@(x) any(strcmp(x,{'point','index'})))
addParameter(p,'flip_xy',false,@(x) islogical(x) || x==1 || x==0)
addParameter(p,'dim',3,@(x) x==2 || x==3)
parse(p,'filename','points',varargin{:})

nPoints = size(points,1);

if p.Results.flip_xy == true
    points(:,1:2) = -points(:,1:2);
end

% Write a file with the points in the transformix file format
fid = fopen(filename,'w');
fprintf(fid,'%s\n',p.Results.type);
fprintf(fid,'%d\n',nPoints);
for k = 1 : nPoints
    if p.Results.dim == 2
        fprintf(fid,'%f %f\n',points(k,:));
    elseif p.Results.dim == 3
        fprintf(fid,'%f %f %f\n',points(k,:));
    end
end
fclose(fid);
end % of the function

