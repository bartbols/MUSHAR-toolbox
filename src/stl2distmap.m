function stl2distmap(stlFile,varargin)
%STL2DISTMAP Calls ShapeWorks and Convert3D to create a signed/absolute
%distance map from an STL file.
%
% Tested to work with ShapeWorks v6.2.1. May not work with earlier
% versions.
%
% Expects ShapeWorks and Convert3D to be installed and added to the path.
%
% INPUT:
% stlFile   : filename of the stl-file from which a distance map will be
%             created.
%
% OPTIONAL INPUT:
% distmap_signed: filename of signed distance map. Default: same name as STL-file but
%                 with extension '_signed.nii.gz'.
% padding       : padding around mesh (determines bounding box of distance
%                 map image). Default: 5 (same units as mesh)
% spacex        : voxel size of distance map in x-direction. Default : 1
% spacey        : voxel size of distance map in y-direction. Default : 1
% spacez        : voxel size of distance map in z-direction. Default : 1

% distmap_abs:    filename of absolute distance map.
% NOTE: if distmap_signed and distmap_abs are both empty, only a signed
%              distance map with the default filename will be created.
%
% OUTPUT: no outputs but file(s) will be created.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% 29/3/2021
% 
% change log:
% 3/2/2022: only allows 'padding' and 'space' as optional arguments to make
% the file compatible with ShapeWorks v6.2.1 ('origin' and and 'size' are no
% longer available as input arguments to mesh-to-dt)

p = inputParser;
addRequired(p,'stlFile',@(x) exist(x,'file')==2)
addParameter(p,'distmap_signed',[])
addParameter(p,'distmap_abs',[])
addParameter(p,'spacex',1,@isscalar) % voxel size of distance map in x-direction
addParameter(p,'spacey',1,@isscalar)
addParameter(p,'spacez',1,@isscalar)
addParameter(p,'padding',5,@isscalar)
parse(p,stlFile,varargin{:});

stlFile        = p.Results.stlFile;
distmap_signed = p.Results.distmap_signed;
distmap_abs    = p.Results.distmap_abs;

delete_signed_distmap = false;
if isempty(distmap_signed)
    % If no filename of a signed distance map is provided, create a signed
    % distance map with the default filename.
    [a,b,c] = fileparts(stlFile);
    distmap_signed = fullfile(a,[b '_signed.nii.gz']);
    if ~isempty(distmap_abs)
        delete_signed_distmap = true;
    end
end

% Call ShapeWorks to build the distance map.

% Build shapeworks command
shapeworks_cmd = sprintf('shapeworks read-mesh --name="%s" mesh-to-dt --sx=%.3f --sy=%.3f --sz=%.3f --pad=%.3f writeimage --name="%s"',...
    stlFile,p.Results.spacex,p.Results.spacey,p.Results.spacez,p.Results.padding,distmap_signed);


if exist(fileparts(distmap_signed),'dir')~=7
    mkdir(fileparts(distmap_signed))
end
[status,cmdout] = system(shapeworks_cmd);

if contains(cmdout,'is not recognized as an internal or external command')
    error('shapeworks not found as an external command. Please install ShapeWorks and try again.')
end

% Flip x and y direction to make distance map compatible with the NIfTI
% coordinate system, and flip the sign for positive values to be outside
% and negative inside the surface.
c3d_cmd = sprintf('c3d "%s" -flip xy -scale -1.0 -o "%s"',...
    distmap_signed,distmap_signed);
[status,cmdout] = system(c3d_cmd);


% Create absolute distance map.
if ~isempty(distmap_abs)
    if exist(fileparts(distmap_abs),'dir')~=7
        mkdir(fileparts(distmap_abs))
    end
    c3d_cmd = sprintf('c3d "%s" -dup -times -sqrt -o "%s"',...
        distmap_signed,distmap_abs);
    [status,cmdout] = system(c3d_cmd);
    if delete_signed_distmap == true
        delete(distmap_signed);
    end
end


end % of function

