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
        N = faceNormal(surf_model);
    end
    
    % Scale the normal vector by half the average edge length.
    edge12 = surf_model.Points(surf_model.ConnectivityList(:,1),:) - surf_model.Points(surf_model.ConnectivityList(:,2),:);
    sf = mean(sqrt(sum(edge12.^2,2)));
    C = squeeze(mean(reshape(surf_model.Points(surf_model.ConnectivityList,:)',3,size(surf_model.ConnectivityList,1),3),3))';
    h_norm = quiver3(C(:,1),C(:,2),C(:,3),...
        sf/2*N(:,1),...
        sf/2*N(:,2),...
        sf/2*N(:,3),0,...
        'Color',p.Results.FaceColor,...
        'AutoScale','off');
end

% Turn on some light, if no light objects exist
if isempty(findobj(gca,'Type','light'))
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

