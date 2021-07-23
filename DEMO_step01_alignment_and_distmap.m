% This demo scripts will rigidly align a collection of muscle shapes and
% then create distance maps from the aligned surfaces. The distance maps
% can subsequently be registered to establish  point-to-point 
% correspondence on the muscle surface and inside the muscle volumes (see
% DEMO_step02).

clear
addpath(genpath('bin'),genpath('src'))
dataFolder = fullfile('example_data');

% 01 = RAND05_VM
% 02 = RAND23_VM 
% 03 = RAND06_VM
% 04 = RAND18_VM
% 05 = RAND10_VM
% 06 = RAND17_VM
%  pairs: 1-2, 3-4, 5-6

%% Rigidly align all surfaces to a reference

% Get the filenames of all STL surfaces in the selected folder.
surfaceFolder = fullfile(dataFolder,'surface');
stl_list      = dir(fullfile(surfaceFolder,'original','*.stl'));
n_surf        = length(stl_list);
fprintf('%d surfaces found in %s\n',n_surf,surfaceFolder)

% Rigidly align all models to a reference surface. The function
% rigid_align_surface uses the iterative closest point algorithms built
% into ShapeWorks.

% Select which number in the list of surfaces is used as the reference.
if exist(fullfile(surfaceFolder,'aligned'),'dir')~=7;mkdir(fullfile(surfaceFolder,'aligned'));end
ref_nr = 1;
refSurface = fullfile(surfaceFolder,'original',stl_list(ref_nr).name); 
fprintf('Surfaces are rigidly aligned to %s.\n',stl_list(ref_nr).name)

for nr = 1 : n_surf
    originalSurface = fullfile(surfaceFolder,'original',stl_list(nr).name); 
    alignedSurface  = fullfile(surfaceFolder,'aligned',stl_list(nr).name);
    transformationFilename = fullfile(surfaceFolder,'aligned',strrep(stl_list(nr).name,'.stl','_transform.mat'));
    if nr == ref_nr
        % The reference model does not need to be aligned. Copy the
        % original model to the aligned folder.
        T = eye(4);
        copyfile(originalSurface,alignedSurface);
        save(transformationFilename,'T')
    else
        % Align model    
        [FV_aligned,T ] = rigid_align_surface(refSurface,originalSurface,alignedSurface);
        % Save the transformation matrix.
        save(transformationFilename,'T')
    end
end
%% Plot surfaces before and after alignment
figure('Color','w')
colors = linspecer(n_surf);
for nr = 1 : n_surf

    subplot(1,2,1);hold on
    plotSurface(fullfile(surfaceFolder,'original',stl_list(nr).name),...
        'FaceColor',colors(nr,:))
    if nr == 1
        title('before alignment')
    end
    subplot(1,2,2);hold on
    plotSurface(fullfile(surfaceFolder,'aligned',stl_list(nr).name),...
        'FaceColor',colors(nr,:))
    if nr == 1
        title('after alignment')
    end
end
hlink = linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim','ZLim','View'});
view(0,15)

%% Make distance maps of all surfaces.

% Make a list of aligned surfaces.
stl_list      = dir(fullfile(surfaceFolder,'aligned','*.stl'));

% Create folders for signed and absolute distance maps.
distmapdir_signed = fullfile(dataFolder,'distmap','signed');
distmapdir_abs    = fullfile(dataFolder,'distmap','abs');
if exist(distmapdir_signed,'dir')~=7;mkdir(distmapdir_signed);end
if exist(distmapdir_abs,'dir')~=7;mkdir(distmapdir_abs);end

% ---------------------------------------------------------------------- %
% Set the bounding box for all distance maps as the maximum bounding box of
% all shapes, plus some padding.
padding = 20; % in mm
voxelsize = [1 1 1]; % voxelsize of distance map in mm

maxvalues = NaN(length(stl_list),3);
minvalues = NaN(length(stl_list),3);
for nr = 1 : length(stl_list)
    % Load each surface and keep the min/max vertex coordinates.
    surfModel = stlread(fullfile(stl_list(nr).folder,stl_list(nr).name));
    maxvalues(nr,:) = max(surfModel.Points);
    minvalues(nr,:) = min(surfModel.Points);
end

% Get over all min/max coordinates and add some padding.
origin = round(min(minvalues)) - padding;
siz    = round(max(maxvalues) - min(minvalues)+2*padding) ./ voxelsize;
save(fullfile(dataFolder,'distmap','dimensions.mat'),'origin','siz');
% ----------------------------------------------------------------------

% Make distance maps for all shapes
% Use the same bounding box for all distance maps.
for nr = 1 : n_surf
    % Build up filenames for distance maps.
    distmap_signed = strrep(stl_list(nr).name,'.stl','_signed.nii.gz');
    distmap_abs    = strrep(stl_list(nr).name,'.stl','_abs.nii.gz');
    fprintf('Making distance map %d of %d (%s)... ',nr,n_surf,distmap_signed)
    surfModel = stlread(fullfile(stl_list(nr).folder,stl_list(nr).name));
    
    % Call ShapeWorks, wrapped in the function stl2distmap, to make the
    % signed and the absolute distance map.
    stl2distmap(...
        fullfile(stl_list(nr).folder,stl_list(nr).name),...
        'distmap_signed',fullfile(distmapdir_signed,distmap_signed),...
        'distmap_abs',fullfile(distmapdir_abs,distmap_abs),...
        'originx', origin(1),...
        'originy', origin(2),...
        'originz', origin(3),...
        'sizex',siz(1),...
        'sizey',siz(2),...
        'sizez',siz(3),...
        'spacex',voxelsize(1),...
        'spacey',voxelsize(2),...
        'spacez',voxelsize(3));
    fprintf(' completed.\n')
end

% Show the middle transverse slice of the last generated signed distance map.
D = load_untouch_nii(fullfile(distmapdir_signed,distmap_signed));
figure
imagesc(D.img(:,:,round(siz(3)/2)),[-1 1]*15)
colormap(redblueTecplot)
colorbar
axis equal off
title('Example slice of distance map')

