% This demo script will register distance maps to a reference shape.
% Registration is done in two steps: in step 1, all distance maps
% are registered to a chosen reference shape and in step 2, all distance
% maps are registered to the mean shape after step 1.
%
% This step requires elastix to be installed and recognised as an external
% command.
% NOTE: DEMO_step01_alignment_and_distmap.m should have been run before
%       running this scripts.
clear
addpath(genpath('bin'),genpath('src'))

dataFolder    = fullfile('example_data');
atlasFolder   = fullfile(dataFolder,'atlas');
distmapFolder = fullfile(dataFolder,'distmap','abs');

%% Registration step 1: Register distance maps to reference scan

% The result of this step is the first estimate of the mean shape and the
% mapping from the chosen reference model to all other models.
refDistmap = fullfile(distmapFolder,'01_abs.nii.gz');
refModel = fullfile(dataFolder,'surface','aligned','01.stl');

% Set parameter files for Elastix: similarity transform followed by a
% b-spline transform with a grid spacing of 10x10x10 mm
parfile = {'parfiles/parameters_Similarity.txt','parfiles/parameters_BSpline10.txt'};

register_distance_maps(...
    distmapFolder,...
    refDistmap,...
    refModel,...
    parfile,...
    'transformFolder',fullfile(dataFolder,'reg2ref','transform'),...
    'resultFolder',fullfile(dataFolder,'reg2ref'));

% Copy the mean shape and make distance maps
if exist(atlasFolder,'dir')~=7
    mkdir(atlasFolder)
end
copyfile(...
    fullfile(dataFolder,'reg2ref','mean_shape.stl'),...
    fullfile(atlasFolder,'atlas.stl'))

% Make distance maps
stl2distmap(...
    fullfile(atlasFolder,'atlas.stl'),...
    'distmap_signed',fullfile(atlasFolder,'atlas_signed.nii.gz'),...
    'distmap_abs', fullfile(atlasFolder,'atlas_abs.nii.gz'),...
    'originx', origin(1),...
    'originy', origin(2),...
    'originz', origin(3),...
    'sizex',siz(1),...
    'sizey',siz(2),...
    'sizez',siz(3),...
    'spacex',voxelsize(1),...
    'spacey',voxelsize(2),...
    'spacez',voxelsize(3));


%% Registration step 2: Register distance maps to the mean shape from step 1
% The result of this step is the mapping from the atlas shape to all other
% shapes
refDistmap    = fullfile(atlasFolder,'atlas_abs.nii.gz');
refModel      = fullfile(atlasFolder,'atlas.stl');

register_distance_maps(...
    distmapFolder,...
    refDistmap,...
    refModel,...
    parfile,...
    'transformFolder',fullfile(atlasFolder,'reg2atlas','transform'),...
    'resultFolder',fullfile(atlasFolder,'reg2atlas'));

%% Visually inspect the result of registration
% Load the registration results
% data.X  : list of corresponding vertex coordinates of all shapes
% data.res_sq_dist: residual error after registration for all vertices of
%                   all shapes
% data.mean_shape : triangulation object with population mean shape
% data.distmap_list : list of filenames of distance maps used for
%                     registration
data = load(fullfile(atlasFolder,'reg2atlas','registration_results.mat'));
n_surf = size(data.X,3);

% NOTE: The transformed reference surface should closely match the original
%       model.
% The residual error is the average distance, across vertices, of the
% transformed reference surface to the original surface (calculated by
% interpolating the absolute distance map of the original surface at the
% vertices of the transformed reference.)

figure('Color','w')
for nr = 1 : n_surf
    subplot(2,3,nr)
    
    % Original (aligned) surface model
    ho = plotSurface(fullfile(dataFolder,'surface','aligned',sprintf('%02d.stl',nr)));
    
    % Reference model transformed to target
    TR = triangulation(data.mean_shape.ConnectivityList, data.X(:,:,nr));
    hr = plotSurface(TR,'FaceColor','none','EdgeColor','k');
    legend([ho hr],{'original','transformed reference'})
    title(sprintf('Residual error: %.2f mm',mean(data.res_sq_dist(:,nr)))) 
end
set(findobj(gcf,'Type','Axes'),'Clipping','off')
axis(findobj(gcf,'Type','Axes'),'off')
hlink = linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim','ZLim','View'});
view(0,15)

%% Show the mean shape and the individual shapes
figure('Color','w');hold on
hm = plotSurface(data.mean_shape,...
    'FaceColor','none',...
    'EdgeColor','k',...
    'LineWidth',2);
view(0,15)
colors = linspecer(n_surf);
h = zeros(1,n_surf);
for nr = 1 : n_surf
    TR = triangulation(data.mean_shape.ConnectivityList,data.X(:,:,nr));
    h(nr) = plotSurface(TR,...
        'EdgeColor',colors(nr,:),'FaceColor',colors(nr,:),...
        'EdgeAlpha',0.5,'FaceAlpha',0.2);
    txt{nr} = sprintf('shape%02d',nr);
end
legend([hm h],[{'mean shape'} txt]);
set(gca,'CLipping','off')
view(0,15)



