function [ points_out ] = transformix_points( points,transform_file,varargin )
%TRANSFORMIX_POINTS transforms 'points' (n x 3 array with physical/voxel
% coordinates) and transforms these points with the Elastix transformation
% file 'transformfile' using Transformix. The transformed points are
% returned.
%
% Bart Bolsterlee, Neuroscience Research Australia
% November 2017
%
% Example:
% points_out = transformix_points(points_in,transform_file)


% Read inputs
p = inputParser;
addRequired(p,'points')
addRequired(p,'transform_file')
addParameter(p,'type','point',@(x) any(strcmp(x,{'point','index'})))
addParameter(p,'flip_xy',true,@(x) islogical(x) || x==1 || x==0)
addParameter(p,'dim',3,@(x) x==2 || x==3)
parse(p,'points','transform_file',varargin{:})

type = p.Results.type;
flip_xy = p.Results.flip_xy;
dim = p.Results.dim;

%%
if nargin < 3
    type = 'point';
end
% Create temporary working directory.
char_list = char(['a':'z' '0':'9']) ;
tmpdir = [];
while exist(tmpdir,'dir') == 7 || isempty(tmpdir)
    tmpdir = fullfile(tempdir,char_list(ceil(length(char_list)*rand(1,8))));
end
mkdir(tmpdir)

try    
        
    % Check if n x 3 or 3 x n array is provided. An array of the same
    % dimensions will be returned.
    
    n = size(points,2);
    if ~ismember(n,[2 3])
        points = points';
    end
    
    if flip_xy == true
        % Write vertices to transformix input points file. Because of 
        % differences in ITK and NIFTI coordinate system, x- and y-values 
        % need to be flipped before transformation.
        if dim == 2
            R = [-1 0;0 -1];
        elseif dim == 3
            R =  [-1 0 0;0 -1 0;0 0 1];
        end
    else
        if dim == 2
            R = eye(2);
        elseif dim == 3 
            R = eye(3);
        end
    end
    
    % Check for NaNs
    nanidx = any(isnan(points),2);   
    
    InputPointsWriter(fullfile(tmpdir,'inputpoints.txt'),...
        points(~nanidx,:)*R,'type',type,'dim',dim)
    
    % Build up the transformix command.
    transformix_cmd = sprintf('transformix -def %s -out %s -tp %s',...
        fullfile(tmpdir,'inputpoints.txt'),...
        tmpdir,...
        transform_file);
    
    % Run transformix
    [status,cmdout] = system(transformix_cmd);
    
    % Read the transformed vertices
    out = OutputPointsReader(fullfile(tmpdir,'outputpoints.txt'),dim);
    
    % Flip vertices back to NIFTI coordinates and save as STL file.
    points_out = NaN(size(points));
    switch type
        case 'point'
            points_out(~nanidx,:) = out.OutputPoint*R;
        case 'index'
            points_out(~nanidx,:) = double(out.OutputIndexFixed)*R;
    end
    if ~ismember(n,[2 3])
        points_out = points_out';
    end
    rmdir(tmpdir,'s')
    
catch ME
    % remove temporary working directory, then throw error message
    rmdir(tmpdir,'s')
    error(getReport(ME))
end

end

