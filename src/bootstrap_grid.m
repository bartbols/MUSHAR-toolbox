function [avg,lb,ub,distr] = bootstrap_grid(X,group1,group2,varargin)
%BOOTSTRAP_ON_GRID calculates the sampling distribution and confidence
%intervals of the mean difference between two groups of values on a grid
%of corresponding nodes.
%
% USAGE:
% [avg,lb,ub,distr] = bootstrap_on_grid(X,group1,group2,varargin)
%
% INPUT
% X   : m x n array of m observations at corresponding nodes in n shapes.
% group 1 : nx1 indices of surfaces in group 1.
% group 2 : mx1 indices of surfaces in group 2.
%
% Optional input arguments, provided as 'argument',<value> pairs:
% N      : number of bootstrapping samples to draw to determine the 95%
%          confidence intervals. Default: N=1000
% paired : if true, the shapes in group 1 and group 2 are paired (e.g.,
%          observations on the same shape before/after some intervention).
% lb_pct : percentile for the lower bound of the confidence interval.
%          Default: 2.5
% ub_pct : percentile for the lower bound of the confidence interval.
%          Default: 97.5
%
% OUTPUT
% avg    : m x 1 array of mean difference per node (group2 minus group1)
% lb     : m x 1 array of lower bound of the confidence interval (default
%          is 2.5th percentile)
% ub     : m x 1 array of lower bound of the confidence interval (default
%          is 97.5th percentile)


p = inputParser;
addRequired(p,'X',@isnumeric)
addRequired(p,'group1',@isnumeric)
addRequired(p,'group2',@isnumeric)
addParameter(p,'N',1000)
addParameter(p,'paired',false)
addParameter(p,'lb_pct',2.5,@isscalar)
addParameter(p,'ub_pct',97.5,@isscalar)
% addParameter(p,'align',false)
parse(p,X,group1,group2,varargin{:});

N = p.Results.N;
% d = NaN(size(X,1)/6,N);
distr = NaN(size(X,1),N);
n1 = numel(group1); % number of shapes in group 1
n2 = numel(group2); % number of shapes in group 2

if p.Results.paired == true
    if n1 ~= n2
        error('Unequal number of samples in groups. For paired-samples bootstrapping the same number of shapes must be in group 1 and 2.')
    end
end
hwait = waitbar(0,'','Name','Bootstrapping grid');
for ii = 1 : N
    waitbar(ii/N,hwait,sprintf('Sample %d of %d',ii,N))
    if p.Results.paired == true % paired samples
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
    distr(:,ii) = nanmean(X(:,sample2) - X(:,sample1),2);    
end

close(hwait)
fprintf('\n')

% Calculate lower and upper bound of the confidence interval.
avg = nanmean(X(:,group2)-X(:,group1),2); % mean difference of actual data
lb = prctile(distr,p.Results.lb_pct,2);   % lower bound of sampling distribution
ub = prctile(distr,p.Results.ub_pct,2);   % upper bound of sampling distribution

end % of function


