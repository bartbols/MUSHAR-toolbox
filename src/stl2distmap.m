function stl2distmap(stlFile,varargin)
%STL2DISTMAP Calls ShapeWorks and Convert3D to create a signed/absolute
%distance map from an STL file.
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
% distmap_abs:    filename of absolute distance map.
% NOTE: if distmap_signed and distmap_abs are both empty, only a signed
%              distance map with the default filename will be created.
% spacex,sizex,originx, etc.  See ShapeWorks documentation for explanation 
% of inputs 'space', 'size' and 'origin'
%
% OUTPUT: no outputs but file(s) be created.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% 29/3/2021

p = inputParser;
addRequired(p,'stlFile',@(x) exist(x,'file')==2)
addParameter(p,'distmap_signed',[])
addParameter(p,'distmap_abs',[])
addParameter(p,'spacex',1)
addParameter(p,'spacey',1)
addParameter(p,'spacez',1)
addParameter(p,'sizex',[])
addParameter(p,'sizey',[])
addParameter(p,'sizez',[])
addParameter(p,'originx',[])
addParameter(p,'originy',[])
addParameter(p,'originz',[])
addParameter(p,'padding',[])
parse(p,stlFile,varargin{:});

stlFile        = p.Results.stlFile;
distmap_signed = p.Results.distmap_signed;
distmap_abs    = p.Results.distmap_abs;
originx        = p.Results.originx;
originy        = p.Results.originy;
originz        = p.Results.originz;
sizex          = p.Results.sizex;
sizey          = p.Results.sizey;
sizez          = p.Results.sizez;
spacex         = p.Results.spacex;
spacey         = p.Results.spacey;
spacez         = p.Results.spacez;

if ~isempty(p.Results.padding)
    fprintf('padding of %.2f is provided. ''origin'' and ''size'' input arguments will be overwritten (if provided)\n',p.Results.padding)
    FV = stlread(stlFile);
    m = min(FV.vertices);
    r = range(FV.vertices);
    pad = p.Results.padding;

    originx = m(1)-pad;
    originy = m(2)-pad;
    originz = m(3)-pad;
    sizex   = ceil((r(1)+2*pad) ./ p.Results.spacex);
    sizey   = ceil((r(2)+2*pad) ./ p.Results.spacey);
    sizez   = ceil((r(3)+2*pad) ./ p.Results.spacez);
end

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

% Build shapeworks command
shapeworks_cmd = sprintf('shapeworks read-mesh --name="%s" mesh-to-dt',...
    stlFile);

% Optional inputs
inputList = {'originx','originy','originz',...
             'sizex','sizey','sizez',...
             'spacex','spacey','spacez'};
for nr = 1 : length(inputList)
    value = eval(inputList{nr});
    if ~isempty(value)     
        shapeworks_cmd = [shapeworks_cmd sprintf(' --%s=%.3f',...
            inputList{nr},value)];    
    end
end
 shapeworks_cmd = [shapeworks_cmd sprintf(' writeimage --name="%s"',...
     distmap_signed)];
 
 % Call ShapeWorks to build the distance map. 
[status,cmdout] = system(shapeworks_cmd);
if contains(cmdout,'is not recognized as an internal or external command')
    error('shapeworks not found as an external commant. Please install ShapeWorks and try again.')
end

% Flip x and y direction to make distance map compatible with the NIfTI
% coordinate system, and flip the sign for positive values to be outside
% and negative inside the surface.
c3d_cmd = sprintf('c3d "%s" -flip xy -scale -1.0 -o "%s"',...
    distmap_signed,distmap_signed);
[status,cmdout] = system(c3d_cmd);


% Create absolute distance map.
if ~isempty(distmap_abs)
    c3d_cmd = sprintf('c3d "%s" -dup -times -sqrt -o "%s"',...
        distmap_signed,distmap_abs);
    [status,cmdout] = system(c3d_cmd);
    if delete_signed_distmap == true
        delete(distmap_signed);
    end
end


end

