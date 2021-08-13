% This demo script will:
% - find corresponding grid points within a group of muscles
% - obtain diffusion tensors at corresponding locations
% - calculate mean fibre orientations for all grid points
% Interpolation and averaging of diffusion tensors is done in the
% log-Euclidean domain. See for more details:
%
% Arsigny, V., Fillard, P., Pennec, X., Ayache, N., 2006. Log-Euclidean
% metrics for fast and simple calculus on diffusion tensors. Magn Reson Med
% 56 (2),p 411-421.http://doi.org/10.1002/mrm.20965

% 01 = RAND05
% 02 = RAND23
% 03 = RAND06
% 04 = RAND18
% 05 = RAND10
% 06 = RAND17

% folder = 'C:\Users\b.bolsterlee\Documents\Working\shape_and_fibre_modelling\TOPHY\data\DTI\RAND05';
% MakeTensor_map(fullfile(folder,'RAND05_DTI_LPCA.fib.gz'),...
%     fullfile(folder,'RAND05_DTI_LPCA.nii.gz'),[0 1],...
%     'example_data/DTI/01_tensor.nii.gz')
clear
dataFolder  = 'example_data';

% Define grid spacing in mm. A regular grid of points will be generated
% inside the reference shape with this spacing. It is set here to a rather 
% coarse grid with 5 mm spacing to have a manageable number for example 
% purposes only. It can be set to lower values for more detailed analyses.
spacing_mm  = 5; 

%% Step 1: make a grid of points in the mean shape
filename.refModel          = fullfile(dataFolder,'registration','step02','ref_shape.stl');
filename.refDistmap_signed = fullfile(dataFolder,'registration','step02','ref_shape_signed.nii.gz');

refModel = stlread(filename.refModel);

% Make a regular grid of points within the reference shape.
[Xg,Yg,Zg] = meshgrid(...
    min(refModel.Points(:,1)):spacing_mm:max(refModel.Points(:,1)),...
    min(refModel.Points(:,2)):spacing_mm:max(refModel.Points(:,2)),...
    min(refModel.Points(:,3)):spacing_mm:max(refModel.Points(:,3)));
gridsize = size(Xg);

% Keep only grid points inside the muscle.
G_ref = [Xg(:) Yg(:) Zg(:)];
dist2surf = interpolate_nii(filename.refDistmap_signed,G_ref);
isInside = dist2surf < 0;
G_ref = G_ref(isInside,:);

% Plot the grid points.
figure('Color','w');hold on
patch('Vertices',refModel.Points,...
    'Faces',refModel.ConnectivityList,...
    'FaceColor',[1 1 1]*0.5,...
    'EdgeColor','none',...
    'FaceAlpha',0.2,...
    'EdgeAlpha',0.5,...
    'Tag','surface',...
    'DiffuseStrength',0.1,...
    'AmbientStrength',0.8,...
    'SpecularStrength',0.1);
plot3(G_ref(:,1),G_ref(:,2),G_ref(:,3),'o',...
    'MarkerSize',6,'MarkerFaceColor','k',...
    'MarkerEdgeColor','y')
axis tight off equal
view(160,45)

light('position',[0 1 0]);light('position',[0 -1 0])
title('Selected points in the reference shape')

%% Get diffusion tensor at corresponding locations
% The diffusion tensor files are expected to be stored in a NIfTI
% file with 4 dimensions in which the first 3 dimensions indicate the size
% of the DTI scan (image dimensions) and the fourth dimension is 6. Each
% component (or stack) in the fourth dimension is one component of the
% diffusion tensor:
% 1 = txx
% 2 = tyy
% 3 = tzz
% 4 = txy
% 5 = txz
% 6 = tyz
% where txx is component (t(1,1) of the full tensor, etc.). Note that only
% 6 values are needed to describe the full tensor because diffusion tensors
% are symmetric (txy = tyx, txz=tzx, tyz=tzx).

tensorFolder = fullfile(dataFolder,'DTI'); 
transformFolder = fullfile(dataFolder,'registration','step02','transform');

ls = dir(fullfile(tensorFolder,'*.nii.gz'));
nShapes = size(ls,1);
nGridPoints = size(G_ref,1);
TENSOR      = NaN([nGridPoints, 6, nShapes]); % matrix with tensor in vector format for all shapes and grid points
G_glob      = NaN([nGridPoints, 3, nShapes]); % matrix with grid locations in global (MRI) coordinate system of each shape
G_aligned   = NaN([nGridPoints, 3, nShapes]); % matrix with grid locations in aligned coordinate system of each shape


for nr = 1 : nShapes
    
    fprintf('Getting tensors for shape %d... ', nr)
    
    % Transform the reference grid to the current muscle using the previously
    % calculated b-spline transform. Note that this will give the grid 
    % locations in the coordinate system of the muscle *after* rigid
    % alignment.
    filename.transform = fullfile(transformFolder,sprintf('ref_shape_abs_to_%02d_abs.txt',nr));
    G_aligned(:,:,nr) = transformix_points(G_ref,filename.transform); 
    
    % Transform to the coordinate *before* rigid alignment, i.e., to the
    % original MRI coordinates.
    filename.rigid_transform = fullfile(dataFolder,'surface','aligned',sprintf('%02d_transform.mat',nr));
    load(filename.rigid_transform); % will load variable 'T', a 4x4 transformation matrix from the aligned to MRI coordinates
    G_glob(:,:,nr) = G_aligned(:,:,nr) * T(1:3,1:3)' + T(1:3,4)';
    
    % Interpolate the tensor map. Interpolation is done in the
    % log-Euclidean domain, so the diffusion tensor first needs to be
    % log-transformed.
    
    % Load tensor file.
    filename.tensor = fullfile(tensorFolder,sprintf('%02d_tensor.nii.gz',nr)); 
    tensor          = load_untouch_nii(filename.tensor); 
    
    % Log transform the tensor (see logTensor.m for details). This is the
    % step that takes the longest to compute.
    tensor.img      = logTensor(vec2tensor(tensor.img)); 
    
    % Interpolate tensor in vector format. Once in log-Euclidean domain,
    % normal interpolation can be used.
    logvec_ip = interpolate_nii(tensor,G_glob(:,:,nr),'linear');
    
    % Transform back to Euclidean domain.
    tensor_ip = expTensor(logvec_ip);
    
    % NOTE: Some programs store diffusion tensors in the ITK coordinate
    % system, which is 180 degr rotated around the z-axis (i.e. flipped x
    % and y-axis) compared to NIfTI coordinate system. If that is the case
    % (as with this example data), the tensor will have to be rotated from 
    % ITK to NIfTI coordinates.
    tensor_ip = rotateTensor(tensor_ip,rotz(pi));
        
    % The tensor is now stored in original MRI coordinates. Rotate to 
    % the reference shape coordinates using the rigid rotation used to
    % align shapes (before nonrigid registration).
    tensor_loc = rotateTensor(tensor_ip,T(1:3,1:3)');
    
    % Store tensor at corresponding locations in vector format.
    TENSOR(:,:,nr) = tensor2vec(tensor_loc);
    fprintf(' completed.\n')
end

% Save the grid locations and the tensor.
save(fullfile(dataFolder,'TENSORDATA.mat'),'TENSOR','G_glob','G_aligned','G_ref','spacing_mm')
fprintf('Diffusion tensors saved as %s.\n',fullfile(dataFolder,'TENSORDATA.mat'))

%% Example: Calculate mean architecture of the group
% First, log-transform the data
logS = NaN(size(TENSOR));
for nr = 1 : size(TENSOR,3) % loop over muscles in the dataset
    logS(:,:,nr) = logTensor(vec2tensor(TENSOR(:,:,nr)));
end

% Once in log-Euclidean domain, the mean can simply be calculated by
% averaging tensors across muscles.
mean_logS = nanmean(logS,3); % mean tensor in log-Euclidean domain (in 1x6 vector format) for all nodes

% Transform back to Euclidean space
meanS = expTensor(mean_logS); % mean tensor in Euclidean domain (in 3x3 tensor format) for all nodes.

%% Example: Plot the individual and group-averaged fibre orientations.
figure('Color','w')
colors = linspecer(size(TENSOR,3));

% Plot individual and group-averaged muscle fibre orientations (primary
% eigenvector of the diffusion tensors)
subplot(1,2,1);hold on
patch('Vertices',refModel.Points,...
    'Faces',refModel.ConnectivityList,...
    'FaceColor','none',...
    'EdgeColor',[1 1 1]*0.5,...
    'EdgeAlpha',0.2);

for nr = 1 : size(TENSOR,3)
    plotTensor(vec2tensor(TENSOR(:,:,nr)),...
        'X',G_ref(:,1),'Y',G_ref(:,2),'Z',G_ref(:,3),...
        'style','stick',...
        'color',colors(nr,:),...
        'scale',spacing_mm*0.8,...
        'ratio',20)
    
end

% Group-averaged fibre orientation
hm = plotTensor(meanS,...
    'X',G_ref(:,1),'Y',G_ref(:,2),'Z',G_ref(:,3),...
    'style','stick',...
    'color',[0 0 0],...
    'scale',spacing_mm,...
    'ratio',8,...
    'facealpha',0.3);
title('Individual (coloured) and group-averaged (black) fibre orientations')

subplot(1,2,2)
plotTensor(meanS,...
    'X',G_ref(:,1),'Y',G_ref(:,2),'Z',G_ref(:,3),...
    'style','stick',...
    'color','ev1',...
    'scale',spacing_mm)
set(findobj(gcf,'Type','Axes'),'Clipping','off')
title('Group-averaged fibre orientations colour-coded to direction')
hlink = linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim','ZLim','View'});
view(140,15)
updateLight

%% Example: Plot the group-averaged diffusion tensors as ellipsoid

figure('Color','w')
hold on
patch('Vertices',refModel.Points,...
    'Faces',refModel.ConnectivityList,...
    'FaceColor','none',...
    'EdgeColor',[1 1 1]*0.5,...
    'EdgeAlpha',0.2);


% Group-averaged diffusion tensors
handle.ellipsoid = plotTensor(meanS,...
    'X',G_ref(:,1),'Y',G_ref(:,2),'Z',G_ref(:,3),...
    'style','ell',...
    'color','ev1',...
    'scale',spacing_mm/4,...
    'facealpha',1.0);
title('Group-averaged diffusion tensors')
