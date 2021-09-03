% This demo script will statistically analyze changes in mean shape between
% groups.

% The example dataset contains 6 vastus medialis muscles from 3 participants
% imaged before and after a period of strength training. This script shows
% an example of how the effect of training on mean muscle shape can be
% analysed.

% shape 1: subject 1, before training
% shape 2: subject 1, after training
% shape 3: subject 2, before training
% shape 4: subject 2, after training
% shape 5: subject 3, before training
% shape 6: subject 3, after training
clear
dataFolder  = fullfile(pwd,'example_data');

group1 = [1 3 5]; % list of shapes before training
group2 = [2 4 6]; % list of shapes after training (order should match group 1 for paired comparisons)

% Load point-correspondence data from registration results
results = load(fullfile(dataFolder,'registration','step02','registration_results.mat'));

% Roll out vertices per shape into a vector so that each column of Xs
% contains all coordinates (x,y and z) of the vertices of a shape. The top
% third of rows contains x-coordinates, the middle third y and the bottom
% third z.
Xs = reshape(results.X,[],size(results.X,3));

%% Bootstrap change in shape.
% The effect of training is defined as the distance from the mean surface
% after training to any point on the surface of the mean surface before
% training (i.e., the point-to-surface distance rather than point-to-point
% distance). The sign of the distance indicates if the point after training
% is located outside the mean surface before training (positive, local 
% expansion) or inside (negative, local contraction). The algorithm used
% to calculate the point-to-surface distance assumes that the normal 
% vectors to the surface point outside, which needs to be verified.

% Verify that the normal vector of the mean shape points outside.
figure('Color','w')
plotSurface(results.mean_shape,'ShowNormals',true)
title('The surface normals should points OUTSIDE')
view(0,15)

% Bootstrap the surface to determine the significance of the effect of
% training on local change in shape. In essence, the bootstrapping approach
% replicates the study N times by randomly sampling (with replacement) from
% the existing shapes. For each bootstrapping replicate, the signed change
% in shape is determined for all vertices. This will give for each vertex
% a distribution of N values for local changes in shape. If the 2.5th
% percentile of this distribution is above 0, the local expansion in shape
% is considered statistically significant. If the 97.5th percentile of the
% distribution is below 0, the local contraction of the shape is considered
% significant. Otherwise, it is not significant (i.e. 95% confidence
% interval includes 0).

% Note that this example only contains 3 pairs of muscles (6 muscles in
% total), which is not sufficient for any statistical analysis. So this is
% just for demonstration purposes only.
bs_results = bootstrap_surface(Xs,...       % 2D array of all vertices (rows) of all shapes (column) - see bootstrap_surface.m for formatting
    results.mean_shape.ConnectivityList,... % list of connectivity (faces) of vertices in Xs
    group1,group2,... % lists of indices of shapes in group 1 and group2 (the effect is group2 minus group1)
    'plot',true,...   % plot mean, lower bound, upper bound and significance map of change in shape
    'paired',true,... % set to true if the data are paired - make sure entries in group 1 and group 2 are corresponding
    'maxdist',5,...   % maxdist is just for plotting purposes
    'N',1000,...       % number of bootstrap replicates. If not provided, the default of 1000 is used.
    'align',true,...  % if true, all shapes are rigiidly aligned with Procrustes alignment before bootstrapping.
    'flipnormals',false,... % set to true if the normal vectors are pointing inside.
    'method','p2s'); % calculate point-to-surface distance
