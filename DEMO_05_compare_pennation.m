% This demo script will statistically analyze changes in pennation angle at
% all nodes in a grid.

% The example dataset contains 6 vastus medialis muscles from 3 participants
% imaged before and after a period of strength training. This script shows
% an example of how the effect of training on local pennation angle can be
% analysed.

% shape 1: subject 1, before training
% shape 2: subject 1, after training
% shape 3: subject 2, before training
% shape 4: subject 2, after training
% shape 5: subject 3, before training
% shape 6: subject 3, after training

clear
dataFolder  = 'example_data';

group1 = [1 3 5]; % list of shapes before training
group2 = [2 4 6]; % list of shapes after training (order should match group 1 for paired comparisons)

% Load diffusion tensors at corresponding locations (created by DEMO_03).
TENSORDATA = load(fullfile(dataFolder,'TENSORDATA'));

exclude = true;
FA_range = [0.1 0.5];

%% Rotate tensor data to its local coordinate system.
% There are different ways to calculate a local coordinate system of a
% muscle. Here, a principal component analysis on the nodes of the muscle
% is used. The 1st principal component is assumed to be the longitudinal
% axis of the muscle. Alternatively, the long axis could be defined as the
% line connecting the origin and insertion of the muscle.

R = pca(TENSORDATA.G_ref); % use pca on all nodes to determine the long axis of the muscle.
R=R(:,[3 2 1]);                % make the z-axis of the local coordinate system the long axis (1st PC)
if R(3,3)<0;R(:,3)=-R(:,3);end % make sure the z-axis points up (just a convention I prefer)
if det(R)<0;R(:,1)=-R(:,1);end

% Set the centroid of all nodes as the origin of the local coordinate
% system.
O = mean(TENSORDATA.G_ref); 
TENSORDATA.G_ref_loc = (TENSORDATA.G_ref-O) * R;

% Rotate tensor data to the local coordinate system.
TENSORDATA.TENSOR = permute(tensor2vec(rotateTensor(vec2tensor(permute(TENSORDATA.TENSOR,[1 3 2])),R)),[1 3 2]);
TENSORDATA.TENSOR_loc = TENSORDATA.TENSOR;

%% Optionally, exclude values based on diffusion properties
if exclude==true
    tensor_in = vec2tensor(permute(TENSORDATA.TENSOR,[1 3 2]));
    [tensor_out,isOutOfRange,bounds] = excludeTensor(tensor_in,...
        'MD',2,...      % exclude tensors with mean diffusivities more than 2 standard deviations away from the mean
        'FA',FA_range); % only include tensors if the FA is between 0.1 and 0.5

    TENSORDATA.TENSOR = permute(tensor2vec(tensor_out),[1 3 2]);
end

%% Calculate and plot mean fibre orientation in group 1 and 2
% Transform tensor into log-Euclidean domain

logS = NaN(size(TENSORDATA.TENSOR));
for nr = 1 : size(TENSORDATA.TENSOR,3) % loop over muscles in the dataset
    logS(:,:,nr) = logTensor(vec2tensor(TENSORDATA.TENSOR(:,:,nr)));
end

% Calculate mean per group, and transform back to Euclidean space.
meanS1 = expTensor(nanmean(logS(:,:,group1),3));
meanS2 = expTensor(nanmean(logS(:,:,group2),3));

% Extract fibre orientation (primary eigenvector) per group
EV1_1 = tensor2ev(meanS1);
EV1_2 = tensor2ev(meanS2);

% Calculate angle in 3D between fibre orientations for all nodes.
mean_ang_diff = acosd(abs(sum(EV1_1 .* EV1_2,2)));
fprintf('Mean (SD across nodes) 3D angular difference between groups is %.2f (%.2f)\n',...
    nanmean(mean_ang_diff),nanstd(mean_ang_diff))

%% Bootstrap for significant differences in pennation angle
% Pennation angle is defined here as the angle of the primary eigenvector
% at a node relative to the long axis (z-axis) of the muscle.

% Calculate the primary eigenvector for all nodes of all muscles.
EV1 = tensor2ev(vec2tensor(permute(TENSORDATA.TENSOR_loc,[1 3 2])));

% Calculate the pennation angle (array with all nodes (rows) for all
% muscles (columns)).
pen_angle = squeeze(acosd(abs(EV1(:,:,3))));

% Bootstrap the pennation angle to determine the significance of the effect
% of training on local pennetion angle. In essence, the bootstrapping
% approach replicates the study N times by randomly sampling (with
% replacement) from the existing muscles. For each bootstrapping replicate,
% the signed change in pennation angle for all nodes is determined. This
% will give for each nodes a distribution of N values for local changes in
% pennation. If the 2.5th percentile of this distribution is above 0, the
% local increase in shape is considered statistically significant. If the
% 97.5th percentile of the distribution is below 0, the local decrease in
% pennation is considered significant. Otherwise, it is not significant
% (i.e. 95th confidence interval includes 0).
%
% Note 1: This example only contains 3 pairs of muscles (6 muscles in
% total), which is not sufficient for any statistical analysis. So this is
% for demonstration purposes only.
%
% Note 2: This bootstrapping approach may not appropriately take into
% account clustering of observations (i.e. the spatial relationship between
% nodes). Other approaches like statistical parametric mapping may be more
% appropriate.

bs_results = bootstrap_grid(pen_angle,...
    group1,group2,...
    'paired',true,...
    'N',1000);

%% Plot results
% Plot group-averaged architectures in the local reference space.

% Average architecture in group 1
figure('Color','w')
subplot(1,4,1);hold on
plotTensor(meanS1,...
    'X',TENSORDATA.G_ref_loc(:,1),...
    'Y',TENSORDATA.G_ref_loc(:,2),...
    'Z',TENSORDATA.G_ref_loc(:,3),...
    'color','ev1',...
    'style','line',...
    'scale',5)
title('Group-averaged architecture group 1')
    
% Average architecture in group 2
subplot(1,4,2);hold on
plotTensor(meanS2,...
    'X',TENSORDATA.G_ref_loc(:,1),...
    'Y',TENSORDATA.G_ref_loc(:,2),...
    'Z',TENSORDATA.G_ref_loc(:,3),...
    'color','ev1',...
    'style','line',...
    'scale',5)
title('Group-averaged architecture group 2')

% Average architecture in group 1 and 2 overlaid.
subplot(1,4,3);hold on
h1 = plotTensor(meanS1,...
    'X',TENSORDATA.G_ref_loc(:,1),...
    'Y',TENSORDATA.G_ref_loc(:,2),...
    'Z',TENSORDATA.G_ref_loc(:,3),...
    'color',[1 0 0],...
    'style','line',...
    'scale',5);
h2 = plotTensor(meanS2,...
    'X',TENSORDATA.G_ref_loc(:,1),...
    'Y',TENSORDATA.G_ref_loc(:,2),...
    'Z',TENSORDATA.G_ref_loc(:,3),...
    'color',[0 0 0],...
    'style','line',...
    'scale',5);
% legend([h1 h2],{'average of group 1','average of group 2'})
title('Group-averaged architecture of group 1 (red) and group 2 (black)')

% Average architecture color-coded to difference in pennation angle.
subplot(1,4,4)
plotTensor(meanS1,...
    'X',TENSORDATA.G_ref_loc(:,1),...
    'Y',TENSORDATA.G_ref_loc(:,2),...
    'Z',TENSORDATA.G_ref_loc(:,3),...
    'color','custom',...
    'cdata',bs_results.avg,...
    'style','line',...
    'scale',5);

colormap(gca,redblueTecplot)
set(gca,'CLim',[-1 1]*10);
title('Difference in pennation angle')
colorbar(gca,'east')

set(findobj(gcf,'Type','Axes'),'Clipping','off')
hlink = linkprop(findobj(gcf,'Type','Axes'),{'XLim','YLim','ZLim','View'});
view(140,15)
updateLight
