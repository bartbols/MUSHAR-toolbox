function results = bootstrap_surface(X,conn,group1,group2,varargin)
%BOOTSTRAP_SURFACE Calculates the mean and confidence interval of the
% distance between two groups of surfaces. The distance is defined as the
% distance of the averaged surface of group 2 to the averaged surface of
% group 1. Bootstrapping is used to calculate the sampling distribution and
% lower/upper bounds.
%
% USAGE:
% results = boostrap_surface(X,conn,group1,group2,varargin)
%
% INPUT
% X   : m x n array of m corresponding points on a selection of n surfaces.
%       Each column m is expected to be the rolled-out vector of a p x 3
%       array of p vertices.
% conn    : ConnectivityList of shape.
% group 1 : nx1 indices of surfaces in group 1.
% group 2 : mx1 indices of surfaces in group 2.
%
% Optional input arguments, provided as 'argument',<value> pairs:
% N      : number of samples from each group to draw to determine the 95%
%          confidence intervals. Default: N=1000
% paired : if true, the shapes in group 1 and group 2 are paired (e.g.,
%          observations on the same shape before/after some intervention).
%          bootstrapping will then sample from the difference vectors
%          between paired shapes. Default: true.
% plot   : if true, a figure with surface coloured by the mean, the lower
%          bound, the upper bound, and a 'significance' (CI crosses 0).
%          Default: true.
% print  : if true, some summary statistics are printed to the screen.
%          Default: true
% maxdist: maximum distance for the colorbar (just for plotting purposes).
% binsize: bin size for histogram. Default: 1
% align  : if true, all shapes are rigidly aligned to the group mean using
%          procrustes alignment
% biascorr : if true, the mean of the sampling distribution is shifted to
%            the mean of the original dataset to reduce sampling bias in
%            the estimates of the confidence interval.
% lb_pct : percentile for the lower bound of the confidence interval.
%          Default: 2.5
% ub_pct : percentile for the lower bound of the confidence interval.
%          Default: 97.5
% proj   : use projection of displacement vector on normal vector.
%          Default: false
% flipnormals: if true, the normals of the surface model are flipped in
%              direction (will only affect the sign of the distance between
%              corresponding points). Default: false.
%
% OUTPUT
% structure with the following fields:
% avg    : m x 1 array with mean effect for all vertices
% lb     : m x 1 array with lower bound for all vertices
% ub     : m x 1 array with upper bound for all vertices
% distr  : m x n array with sampling distribution for all vertices
% ...and some options as described under INPUT
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% July 2021


p = inputParser;
addRequired(p,'X',@isnumeric)
addRequired(p,'conn',@isnumeric)
addRequired(p,'group1',@isnumeric)
addRequired(p,'group2',@isnumeric)
addParameter(p,'N',1000)
addParameter(p,'paired',false)
addParameter(p,'biascorr',true)
addParameter(p,'align',false)
addParameter(p,'flipnormals',false)
addParameter(p,'plot',true)
addParameter(p,'print',true)
addParameter(p,'maxdist',[])
addParameter(p,'binsize',1)
addParameter(p,'lb_pct',2.5,@isscalar)
addParameter(p,'ub_pct',97.5,@isscalar)
addParameter(p,'proj',false)
parse(p,X,conn,group1,group2,varargin{:});


results.options.N           = p.Results.N;   % number of bootstrap replicates
results.options.biascorr    = p.Results.biascorr;
results.options.align       = p.Results.align;
results.options.flipnormals = p.Results.flipnormals;
results.options.lb_pct      = p.Results.lb_pct;
results.options.ub_pct      = p.Results.ub_pct;
results.options.proj        = p.Results.proj;
results.options.paired      = p.Results.paired;

n1 = numel(group1); % number of shapes in group 1
n2 = numel(group2); % number of shapes in group 2
    
if results.options.align == true
    % Pre-align all shapes to the mean shape using procrustes aligment.
    X = rigid_align_to_mean(X);
end

if p.Results.paired == true
    if n1 ~= n2
        error('Unequal number of samples in groups. For paired-samples bootstrapping the same number of shapes must be in group 1 and 2.')
    end
end

if p.Results.flipnormals == true
    % Change order of connectivity list to make surface normals point in
    % the opposite direction.
    conn = conn(:,[2 1 3]);
end


% Calculate shape change of the original data
X1  = reshape(mean(X(:,group1),2),[],3);
X2  = reshape(mean(X(:,group2),2),[],3);
TR1 = triangulation(conn,X1); % mean of group 1
TR2 = triangulation(conn,X2); % mean of group 2
d_vec = TR2.Points - TR1.Points;
N1 = vertexNormal(TR1);
sgn = sign(sum(d_vec .* N1,2));

if p.Results.proj == true
    % Magnitude and direction of projection of displacement vector on normal vector
    projection = d_vec .* N1; % projection of displacement vector on normal vector
    results.avg = sgn .* sqrt(sum(projection.^2,2));
else
    % Magnitude and direction of displacement vector
    results.avg = sgn .* sqrt(sum(d_vec.^2,2));
end

hwait = waitbar(0,'','Name','Bootstrapping surface');
results.distr = NaN(size(X,1)/3,results.options.N); % array with distances for all vertices and samples
for ii = 1 : results.options.N
    waitbar(ii/results.options.N,hwait,sprintf('%d of %d',ii,results.options.N))
    if results.options.paired == true % paired samples
        % Draw a random sample of indices from 1 to n1 (number of shapes
        % per group)
        random_sample = datasample(1:n1,n1);
        
        % Get the indices for paired samples.
        sample1 = group1(random_sample);
        sample2 = group2(random_sample);
        
    else % unpaired samples
        % Sample shapes from each group independently (with replacement)
        sample1 = group1(datasample(1:n1,n1));
        sample2 = group2(datasample(1:n2,n2));
    end
    
    % Calculate the mean shapes in each sampled group.
    X1_sampled = reshape(mean(X(:,sample1),2),[],3);
    X2_sampled = reshape(mean(X(:,sample2),2),[],3);
    
    % Make reference shape and distance map of reference shape.
    TR1 = triangulation(conn,X1_sampled);
    TR2 = triangulation(conn,X2_sampled);
    d_vec = TR2.Points - TR1.Points;
    N1 = vertexNormal(TR1);
    sgn = sign(sum(d_vec .* N1,2));
    if results.options.proj == true
        % Magnitude and direction of projection of displacement vector on normal vector
        projection = d_vec .* N1;
        results.distr(:,ii) = sgn .* sqrt(sum(projection.^2,2));
    else
        % Magnitude and direction of displacement vector
        results.distr(:,ii) = sgn .* sqrt(sum(d_vec.^2,2));
    end
end
close(hwait)
fprintf('\n')

if results.options.biascorr == true
    % Bias correction: correct sampling distribution to have a mean equal to
    % the true mean
    results.bias = nanmean(results.distr,2) - results.avg;
    results.distr = results.distr - results.bias;
    results.avg = results.avg;
end
% Calculate lower and upper bound of the 95% confidence interval. Correct
% for small interpolation errors by subtracting the distance (d0) of the mean
% shape to the distance map. d0 should be zero if the distance map is
% perfect.
results.lb = prctile(results.distr,results.options.lb_pct,2); % lower bound
results.ub = prctile(results.distr,results.options.ub_pct,2); % upper bound

if p.Results.print == true
    % Print results to the screen
    fprintf(' -------- BOOTSTRAP RESULTS --------\n')
    fprintf('Average expansion (SD)                   : %.1f (%.1f)\n',mean(results.avg),std(results.avg))
    fprintf('Significant contraction (%% of vertices) : %.2f \n',sum(results.ub<0)/numel(results.lb)*100)
    fprintf('Not significant         (%% of vertices) : %.2f \n',sum(results.lb<0 & results.ub>0)/numel(results.lb)*100)
    fprintf('Significant expansion   (%% of vertices) : %.2f \n',sum(results.lb>0)/numel(results.lb)*100)
end


%% Plot the results
if p.Results.plot == true
    if isempty(p.Results.maxdist)
        maxdist = ceil(max(abs(results.distr(:))));
    else
        maxdist = p.Results.maxdist;
    end
    
    figure('Color','w')
    hs = zeros(1,4);
    for i = 1 : 4
        hs(i) = subplot(2,4,i);hold on
        hp = patch('Vertices',X1_sampled,...
            'Faces',p.Results.conn,...
            'FaceColor','interp',...
            'EdgeColor','none',...
            'FaceAlpha',1,...
            'EdgeAlpha',0.5,...
            'DiffuseStrength',0.1,...
            'AmbientStrength',0.8,...
            'SpecularStrength',0.1);
        switch i
            case 1 % mean difference
                cdata = results.avg;
                set(hp,'FaceVertexCData',cdata)
                titleTxt = 'mean surface-to-surface distance';
                cmap = redblueTecplot(100);
                set(gca,'CLim',[-1 1]*maxdist)
            case 2 % lower boundary (2.5 percentile)
                cdata = results.lb;
                set(hp,'FaceVertexCData',cdata)
                titleTxt = 'lower bound';
                cmap = redblueTecplot(100);
                set(gca,'CLim',[-1 1]*maxdist)
            case 3 % upper boundary (97.5 percentile)
                cdata = results.ub;
                set(hp,'FaceVertexCData',cdata)
                titleTxt = 'upper bound';
                cmap = redblueTecplot(100);
                set(gca,'CLim',[-1 1]*maxdist)
            case 4 % CI includes 0, >0 or <0
                cdata = ones(size(results.lb))*1;
                cdata(results.ub < 0) = 0; % upper bound < 0 - significant negative change
                cdata(results.lb > 0) = 2; % lower bound > 0 - significant positive change
                %                 cdata = double(lb <= 0 & ub >= 0);
                set(hp,'FaceVertexCData',cdata)
                %                 titleTxt = 'yellow = significant (CI does not include 0)';
                titleTxt = 'significance map';
                cmap = redblueTecplot(3);
                %                 cmap = [1 1 0;0 0 0;1 0 0];
                set(gca,'CLim',[0 2])
        end
        colormap(gca,cmap)
        
        hc = colorbar;
        if i==4
            set(hc,'Ticks',[0 1 2],'TickLabels',{'d<0','NS','d>0'})
        end
        axis tight off equal
        view(-80,25)
        light('position',[0 1 0]);light('position',[0 -1 0])
        
        xlim = get(gca,'Xlim');
        ylim = get(gca,'Ylim');
        zlim = get(gca,'Zlim');
        padding = 5;
        axis(gca,[xlim+[-1 1]*padding ylim+[-1 1]*padding zlim+[-1 1]*padding])
        title(titleTxt)
        
        subplot(2,4,i+4);hold on
        bins = floor(min(cdata)):p.Results.binsize:ceil(max(cdata));
        if isscalar(bins);bins=[-1 0 1]*p.Results.binsize + mean(cdata);end
        if i == 4
            bins = [0 1 2];
        end
        colored_hist(cdata,bins)
        colormap(gca,cmap)
        if i == 4
            set(gca,'XTick',[0 1 2],...
                'XTickLabels',{'sign. contr','N.S.','sign. exp.'},...
                'XLim',[-0.5 2.5],...
                'CLim',[0 2])
        else
            xlabel('Distance (mm)')
            ylabel('Fraction of vertices')
        end
    end
    set(hs,'Clipping','off')
    global hlink
    hlink = linkprop(hs,{'XLim','YLim','ZLim','View'});
end

end % of function


