function register_distance_maps(distmapFolder,refDistmap,refModel,parFile,varargin)
%REGISTER_DISTANCE_MAPS registers distance maps using Elastix.

%% Check input arguments
p = inputParser;
addRequired(p,'distmapFolder',@(x) exist(x,'dir')==7)
% addRequired(p,'stlFolder',@(x) exist(x,'dir')==7)
addRequired(p,'refDistmap',@(x) exist(x,'file')==2)
addRequired(p,'refModel',@(x) exist(x,'file')==2)
addRequired(p,'parFile')
addParameter(p,'transformFolder',[])
addParameter(p,'resultFolder',[])
addParameter(p,'selection',[])
addParameter(p,'overwrite',true)
% addParameter(p,'makeatlas',true,@(x) x==0 || x==1 || islogical(x) )

parse(p,distmapFolder,refDistmap,refModel,parFile,varargin{:});

distmapFolder = p.Results.distmapFolder; % folder with distance maps
% stlFolder     = p.Results.stlFolder;   % folder with stl surfaces
refDistmap    = p.Results.refDistmap;   % reference distance map
refModel      = p.Results.refModel;     % reference surface model
resultFolder  = p.Results.resultFolder; % folder where mean shape and distance maps will be saved to
transformFolder = p.Results.transformFolder;
% makeatlas   = p.Results.makeatlas;   % flag to make an atlas model
parFile     = p.Results.parFile;       % elastix parameter file
selection      = p.Results.selection;  % selection of shapes to register
overwrite      = p.Results.overwrite;  % overwrite transforms if they already exist

%% Register to reference model
if isempty(transformFolder)
    rootFolder      = fileparts(distmapFolder);
    transformFolder = fullfile(rootFolder,'transform');
end

if exist(transformFolder,'dir') ~= 7
    mkdir(transformFolder)
end

% Make list of all nifti files in the distance map folder
distmap_list = dir(fullfile(distmapFolder,'*.nii.gz'));
if ~isempty(selection)
    distmap_list = distmap_list(selection);
end

% Load the reference model
[~,tmp,~]  = fileparts(refDistmap);
[~,short_ref_name] = fileparts(tmp);
ref = stlread(refModel);

% Create NaN arrays.
% Array with all transformed coordinates
X = NaN([size(ref.Points) length(distmap_list)]); 

% Array with squared distance of transformed surface to original surface.
res_sq_dist  = NaN([size(ref.Points,1) length(distmap_list)]); 

% Open a wait bar
h_wait = waitbar(0,'','Name','Registration progress');
for distmap_nr = 1 : length(distmap_list)
    % Update waitbar
    waitbar((distmap_nr-1)/length(distmap_list),...
        h_wait,...
        sprintf('Registering %d of %d',distmap_nr,length(distmap_list)))
    
    % Get filename of 'moving' distance map
    movingDistmap = fullfile(distmapFolder,distmap_list(distmap_nr).name); % distance map of moving model.    
    [~,b,~]  = fileparts(movingDistmap);
    [~,short_moving_name] = fileparts(b);
    
    % The reference model doesn't need to be registered on itself.
    if strcmp(short_ref_name,short_moving_name)
        X(:,:,distmap_nr) = ref.Points;
        continue
    end
    
    % Build up transformation filename
    transform_file = fullfile(transformFolder,sprintf('%s_to_%s.txt',short_ref_name,short_moving_name));
    
    % Register the distance maps using Elastix. Save only the transform
    % file.
    if overwrite == true || exist(transform_file,'file')~=2
        reg_elastix(refDistmap,movingDistmap,parFile,'transform_file',transform_file)
    end
    
    % Calculate the residual distance on the squared distance map of the moving
    % model at the vertices of the transformed reference model. This
    % should be close to 0.
    D_m = load_untouch_nii(movingDistmap);
    X(:,:,distmap_nr) = transformix_points(ref.Points,transform_file);    
    res_sq_dist(:,distmap_nr) = interpolate_nii(D_m,X(:,:,distmap_nr));
    
    % Plot moving model and transformed reference model for a visual
    % check of alignment.
    %         moving = stlread(fullfile(modelPath,[short_moving_name '.stl']));
    %         figure
    %         hold on
    %         axis equal
    %         patch('Vertices',moving.Points,...
    %             'Faces',moving.ConnectivityList,...
    %             'FaceColor','none',...
    %             'EdgeColor','r')
    %
    %         patch('Vertices',X(:,:,distmap_nr),...
    %             'Faces',ref.ConnectivityList,...
    %             'FaceColor','none',...
    %             'EdgeColor','k')
end
close(h_wait) % close the waitbar

% Report results
fprintf('--- Residual distance of transformed reference model to moving model ---\n')
fprintf('%-20s: %s (%s)\n',...
    'Model','mean','SD')
for distmap_nr = 1 : length(distmap_list)
    movingDistmap = fullfile(distmapFolder,distmap_list(distmap_nr).name); % distance map of moving model.
    [~,b,~]  = fileparts(movingDistmap);
    [~,short_moving_name] = fileparts(b);
    fprintf('%-20s: %.2f ( %.2f)\n',...
        short_moving_name,...
        nanmean(sqrt(res_sq_dist(:,distmap_nr))),...
        nanstd(sqrt(res_sq_dist(:,distmap_nr))))
end

%% Write mean shape
if isempty(resultFolder)
    resultFolder = fullfile(fileparts(distmapFolder),'mean');
end

if exist(resultFolder,'dir') ~= 7
    mkdir(resultFolder)
end

% Make atlas model (mean of all registered models)
mean_shape = triangulation(ref.ConnectivityList,nanmean(X,3));

% Write mean shape as an STL file
stlwrite(mean_shape,fullfile(resultFolder,'mean_shape.stl'));
fprintf('Atlas model saved as %s.\n',fullfile(resultFolder,'mean_shape.stl'))

% Save mean models
save(fullfile(resultFolder,'registration_results.mat'),'mean_shape','X','res_sq_dist','distmap_list')


end % of function