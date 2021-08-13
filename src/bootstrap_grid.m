function results = bootstrap_grid(X,group1,group2,varargin)
%BOOTSTRAP_GRID calculates the sampling distribution and confidence
%intervals of the mean difference between two groups of values on a grid
%of corresponding nodes.
%
% USAGE:
% results = bootstrap_grid(X,group1,group2,varargin)
%
% INPUT
% X   : m x n array of m observations at corresponding nodes in n shapes.
% group 1 : nx1 indices of surfaces in group 1.
% group 2 : mx1 indices of surfaces in group 2.
%
% Optional input arguments, provided as 'argument',<value> pairs:
% N      : number of bootstrap replicates to determine the sampling
%          distribution. Default: N=1000
% paired : if true, the shapes in group 1 and group 2 are paired (e.g.,
%          observations on the same shape before/after some intervention).
% lb_pct : percentile for the lower bound of the confidence interval.
%          Default: 2.5
% ub_pct : percentile for the lower bound of the confidence interval.
%          Default: 97.5
% biascorr : if true, the mean of the sampling distribution is shifted to
%            the mean of the original dataset to reduce sampling bias in
%            the estimates of the confidence interval.
% print  : if true, some summary statistics are printed to the screen.
%          Default: true

%
% OUTPUT
% structure with the following fields:
% avg    : m x 1 array of mean difference per node (group2 minus group1)
% lb     : m x 1 array of lower bound of the confidence interval (default
%          is 2.5th percentile)
% ub     : m x 1 array of lower bound of the confidence interval (default
%          is 97.5th percentile)
% distr  : m x n array with sampling distribution for all nodes
% ...and some options as described under INPUT
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% August 2021


p = inputParser;
addRequired(p,'X',@isnumeric)
addRequired(p,'group1',@isnumeric)
addRequired(p,'group2',@isnumeric)
addParameter(p,'N',1000)
addParameter(p,'paired',false)
addParameter(p,'biascorr',true)
addParameter(p,'print',true)
addParameter(p,'lb_pct',2.5,@isscalar)
addParameter(p,'ub_pct',97.5,@isscalar)
parse(p,X,group1,group2,varargin{:});

% Store some options in the results structure.
results.options.N           = p.Results.N;   % number of bootstrap replicates
results.options.biascorr    = p.Results.biascorr;
results.options.lb_pct      = p.Results.lb_pct;
results.options.ub_pct      = p.Results.ub_pct;
results.options.paired      = p.Results.paired;

N = results.options.N;
% d = NaN(size(X,1)/6,N);
results.distr = NaN(size(X,1),N);
n1 = numel(group1); % number of shapes in group 1
n2 = numel(group2); % number of shapes in group 2

if results.options.paired == true
    if n1 ~= n2
        error('Unequal number of samples in groups. For paired-samples bootstrapping the same number of shapes must be in group 1 and 2.')
    end
end
hwait = waitbar(0,'','Name','Bootstrapping grid');
for ii = 1 : N
    waitbar(ii/N,hwait,sprintf('Sample %d of %d',ii,N))
    if results.options.paired == true % paired samples
        % Draw a random sample (with replacement)
        random_sample = datasample(1:n1,n1);
        
        % Get indices of sampled grids
        sample1 = group1(random_sample);
        sample2 = group2(random_sample);
        
    else % unpaired samples
        % Sample indices from each group independently (with replacement)
        sample1 = group1(datasample(1:n1,n1));
        sample2 = group2(datasample(1:n2,n2));
    end 
    % Calculate difference in group means.
    results.distr(:,ii) = nanmean(X(:,sample2) - X(:,sample1),2);    
end

close(hwait)
fprintf('\n')

% Calculate lower and upper bound of the confidence interval.
results.avg = nanmean(X(:,group2)-X(:,group1),2); % mean difference of actual data

if results.options.biascorr == true
    % Bias correction: correct sampling distribution to have a mean equal to
    % the true mean
    results.bias = nanmean(results.distr,2) - results.avg;
    results.distr = results.distr - results.bias;
end

results.lb = prctile(results.distr,results.options.lb_pct,2);   % lower bound of sampling distribution
results.ub = prctile(results.distr,results.options.ub_pct,2);   % upper bound of sampling distribution

if p.Results.print == true
    % Print results to the screen
    fprintf(' -------- BOOTSTRAP RESULTS --------\n')
    % Remove NaNs
    notnan_idx = ~isnan(results.avg);
    n = sum(notnan_idx);
    fprintf('%d nodes were NaN (%.2f%% of total nodes)\n',sum(~notnan_idx),sum(~notnan_idx)/numel(results.avg)*100);
    fprintf('n                                   : %d\n',sum(notnan_idx))
    fprintf('Mean (SD across nodes) group 1      : %.1f (%.1f)\n',mean(nanmean(X(notnan_idx,group1),2)),std(nanmean(X(notnan_idx,group1),2)))
    fprintf('Mean (SD across nodes) group 2      : %.1f (%.1f)\n',mean(nanmean(X(notnan_idx,group2),2)),std(nanmean(X(notnan_idx,group2),2)))
    fprintf('Average difference (SD)             : %.1f (%.1f)\n',nanmean(results.avg(notnan_idx)),nanstd(results.avg(notnan_idx)))
    fprintf('Significantly smaller (%% of nodes)  : %.1f \n',nansum(results.ub(notnan_idx)<0)/n*100)
    fprintf('Not significant       (%% of nodes)  : %.1f \n',nansum(results.lb(notnan_idx)<0 & results.ub(notnan_idx)>0)/n*100)
    fprintf('Significantly larger  (%% of nodes)  : %.1f \n',nansum(results.lb(notnan_idx)>0)/n*100)
end


end % of function


