function varargout = plotSurface( surf_model,varargin )
%PLOTSURFACE Plots a surface model with the default settings, or custom
%settings when provided.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% March 2018
% 
% Usage:
% handle = plotSurface(surf_model,'parameter',<value>)
% 
% surf_model : structure containing fields vertices and faces or filename 
%              of an STL-file.
% Optional inputs:
% FaceColor  : color of the faces of the surface model. Default: 'r'
% EdgeColor  : color of the edge of the surface model. Default: 'none'
% FaceAlpha  : transparency value of the faces. Default: 0.3
% EdgeAlpha  : transparency value of the edges. Default: 1
% LineWidth  : linewidth of the edges. Default: 1
% ShowNormals : Show normal vectors of the surface model.

if nargin == 0
    [file,path] = uigetfile('*.stl','Select an STL-file');
    if file == 0;return;end
    surf_model = fullfile(path,file);
    
end
p = inputParser;
addRequired(p,'surf_model')
addParameter(p,'FaceColor','r')
addParameter(p,'FaceAlpha',0.3,@isscalar)
addParameter(p,'EdgeColor','none')
addParameter(p,'EdgeAlpha',1,@isscalar)
addParameter(p,'LineWidth',1,@isscalar)
addParameter(p,'DiffuseStrength',0.1,@isscalar)
addParameter(p,'AmbientStrength',0.8,@isscalar)
addParameter(p,'SpecularStrength',0.1,@isscalar)
addParameter(p,'ShowNormals',false)

parse(p,surf_model,varargin{:})

% Read the STL-file 
if ischar(surf_model)
    surf_model = stlread(surf_model);
end
    
holdstate = ishold(gca);
hold on
h = patch('Vertices',surf_model.Points,...
               'Faces',surf_model.ConnectivityList,...
               'FaceColor',p.Results.FaceColor,...
               'FaceAlpha',p.Results.FaceAlpha,...
               'EdgeColor',p.Results.EdgeColor,...
               'EdgeAlpha',p.Results.EdgeAlpha,...
               'LineWidth',p.Results.LineWidth,...
               'DiffuseStrength',p.Results.DiffuseStrength,...
               'AmbientStrength',p.Results.AmbientStrength,...
               'SpecularStrength',p.Results.SpecularStrength);

if p.Results.ShowNormals == true
    if ~isfield(surf_model,'normals')
        surf_model.normals = facenormals(surf_model);
    end
    
    C = squeeze(mean(reshape(surf_model.vertices(surf_model.faces,:)',3,size(surf_model.faces,1),3),3))';
    h_norm = quiver3(C(:,1),C(:,2),C(:,3),...
        surf_model.normals(:,1),...
        surf_model.normals(:,2),...
        surf_model.normals(:,3),0,...
        'Color','r',...
        'AutoScale','off');
end

% Turn on some light, if no light objects exist
if isempty(findobj(gcf,'Type','light'))
    light('position',[0 1 0]);
    light('position',[0 -1 0]);
end
axis equal
if holdstate == 0
    hold(gca,'off')
else
    hold(gca,'on')
end

if nargout > 0
    varargout{1} = h;
    if nargout > 1
        varargout{2} = h_norm;
    end
end
end

