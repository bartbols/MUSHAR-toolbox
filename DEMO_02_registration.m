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
dataFolder           = fullfile('example_data');

%% Registration step 1: Register distance maps to reference scan

% The result of this step is the first estimate of the mean shape and the
% mapping from the chosen reference model to all other models.
filename.refDistmap = fullfile(dataFolder,'distmap','abs','01_abs.nii.gz');
filename.refModel   = fullfile(dataFolder,'surface','aligned','01.stl');

% Set parameter files for Elastix: similarity transform followed by a
% b-spline transform with a grid spacing of 10x10x10 mm
% parfile = {'parfiles/parameters_BSpline10.txt'};
parfile = {'parfiles/parameters_BSpline10_tmp.txt'};

results01 = register_distance_maps(...
    fullfile(dataFolder,'distmap','abs'),...
    filename.refDistmap,...
    filename.refModel,...
    parfile,...
    'transformFolder',fullfile(dataFolder,'registration','step01','transform'),...
    'resultFolder',fullfile(dataFolder,'registration','step01'));

%% Step 2: registration to the mean shape of step 1
% Make the mean shape from step 1 the reference shape for step 2 of
% registration. Make distance maps of the mean shape of step 1.

if exist(fullfile(dataFolder,'registration','step02'),'dir')~=7
    mkdir(fullfile(dataFolder,'registration','step02'))
end

filename.refModel          = fullfile(dataFolder,'registration','step02','ref_shape.stl');
filename.refDistmap_abs    = fullfile(dataFolder,'registration','step02','ref_shape_abs.nii.gz');
filename.refDistmap_signed = fullfile(dataFolder,'registration','step02','ref_shape_signed.nii.gz');
copyfile(fullfile(dataFolder,'registration','step01','mean_shape.stl'),...
    filename.refModel)

% Use the same spatial dimensions as the dimension of all the original
% distance maps.
dims = load(fullfile(dataFolder,'distmap','dimensions.mat'));

stl2distmap(...
    refModel,...
    'distmap_signed',filename.refDistmap_signed,...
    'distmap_abs',filename.refDistmap_abs,...
    'originx', dims.origin(1),...
    'originy', dims.origin(2),...
    'originz', dims.origin(3),...
    'sizex',dims.siz(1),...
    'sizey',dims.siz(2),...
    'sizez',dims.siz(3),...
    'spacex',dims.voxelsize(1),...
    'spacey',dims.voxelsize(2),...
    'spacez',dims.voxelsize(3));

% Register distance maps to the mean shape of step 1
results02 = register_distance_maps(...
    fullfile(dataFolder,'distmap','abs'),...
    filename.refDistmap_abs,...
    filename.refModel,...
    parfile,...
    'transformFolder',fullfile(dataFolder,'registration','step02','transform'),...
    'resultFolder',fullfile(dataFolder,'registration','step02'));

%% Visually inspect the result of registration
% Load the registration results
n_surf = size(results.X,3);

% NOTE: The transformed reference surface should closely match the original
%       model.
% The residual error is the average distance, across vertices, of the
% transformed reference surface to the original surface (calculated by
% interpolating the absolute distance map of the original surface at the
% vertices of the transformed reference.)

figure('Color','w','Name','Original and transformed reference')
for nr = 1 : n_surf
    subplot(1,6,nr)
    
    % Original (aligned) surface model
    ho = plotSurface(fullfile(dataFolder,'surface','aligned',sprintf('%02d.stl',nr)),...
        'EdgeColor','r','EdgeAlpha',0.2);
    
    % Reference model transformed to target
    TR = triangulation(results.mean_shape.ConnectivityList, results02.X(:,:,nr));
    hr = plotSurface(TR,'FaceColor','none','EdgeColor','k');
    legend([ho hr],{'original','transformed reference'})
    title(sprintf('Residual error: %.2f mm',mean(results02.res_dist(:,nr)))) 
end
% Set some axes properties
set(findobj(gcf,'Type','Axes'),'Clipping','off')
axis(findobj(gcf,'Type','Axes'),'off')
% Link the viewing angle and axis limits of all subplots.
hlink = linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim','ZLim','View'});
view(0,15)

%% Show the mean shape and the individual shapes
figure('Color','w');hold on
hm = plotSurface(results.mean_shape,...
    'FaceColor','none',...
    'EdgeColor','k',...
    'LineWidth',2);
view(0,15)
colors = linspecer(n_surf);
h = zeros(1,n_surf);
for nr = 1 : n_surf
    TR = triangulation(results.mean_shape.ConnectivityList,results.X(:,:,nr));
    h(nr) = plotSurface(TR,...
        'EdgeColor',colors(nr,:),'FaceColor',colors(nr,:),...
        'EdgeAlpha',0.5,'FaceAlpha',0.2);
    txt{nr} = sprintf('shape%02d',nr);
end
legend([hm h],[{'mean shape'} txt]);
set(gca,'CLipping','off')
view(0,15)



