function [tensor_out,isOutOfRange,bounds] = excludeTensor(tensor,varargin)
%EXCLUDETENSOR Sets values in diffusion tensor field S to NaN if diffusion
%properties are outside the specified boundaries.
%
% [tensor_out,isOutOfRange ] = excludeTensor(tensor,varargin)
% 
% REQUIRED INPUT:
% tensor   : tensor field where the last two dimensions are 3x3
%
% OPTIONAL INPUT: provided as 'property',<value>-pairs:
% MD       : scalar or 1x2 vector with range for mean diffusivity.
%           If a scalar is provided, it is the number of standard 
%           deviations a value is removed from the mean to be excluded.
%           If a 1x2 vector is provided, it is the lower/upper bound of the
%           range.
% L1       : same format as MD, but for exclusion based on primary
%            eigenvalue
% L2       : same format as MD, but for exclusion based on secondar
%            eigenvalue
% L3       : same format as MD, but for exclusion based on tertiary
%            eigenvalue
% L3       : same format as FA, but for exclusion based on fractional
%            anisotropy
% print    : true/false. If true, the number and percentage of values that
%            were excluded (set to NaN) are printed to the screen.
%
% OUTPUT:
% tensor_out : tensor with values outside the provided range for any of the
%              diffusion properties set to NaN.
% isOutOfRange : logical array with, per element of the tensor field, 5
%                values indicating whether the L1, L2 L3, MD and FA were
%                outside the provided range (0/1=inside/outside range).
% bounds       : boundaries used for each of the diffusion properties (L1,
%                L2, L3, MD, FA)

%
% Bart Bolsterlee
% Neuroscience Research Australia
% July 2021
% 

p = inputParser;
addRequired(p,'tensor')
addParameter(p,'MD',[0 Inf])
addParameter(p,'L1',[0 Inf])
addParameter(p,'L2',[0 Inf])
addParameter(p,'L3',[0 Inf])
addParameter(p,'FA',[0 1])
addParameter(p,'print',true)
parse(p,tensor,varargin{:});

% Get input dimensions.
nd = ndims(tensor);
sz = size(tensor);

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-2));
tensor = reshape(tensor,[n,3,3]);

% Get diffusion properties
[L1,L2,L3,MD,FA] = tensor2diffprop(tensor);

% Define thresholds
bounds = NaN(5,2);
if isscalar(p.Results.L1)
        % Input values are min/max number of standard deviations away from
        % the mean
        bounds(1,:) = [-1 1] * p.Results.L1 * nanstd(L1(:)) + nanmean(L1(:));
else
    bounds(1,:) = p.Results.L1;
end
if isscalar(p.Results.L2)
        % Input values are min/max number of standard deviations away from
        % the mean
        bounds(2,:) = [-1 1] * p.Results.L2 * nanstd(L2(:)) + nanmean(L2(:));
else
    bounds(2,:) = p.Results.L2;
end
if isscalar(p.Results.L3)
        % Input values are min/max number of standard deviations away from
        % the mean
        bounds(3,:) = [-1 1] * p.Results.L3 * nanstd(L3(:)) + nanmean(L3(:));
else
    bounds(3,:) = p.Results.L3;
end
if isscalar(p.Results.MD)
        % Input values are min/max number of standard deviations away from
        % the mean
        bounds(4,:) = [-1 1] * p.Results.MD * nanstd(MD(:)) + nanmean(MD(:));
else
    bounds(4,:) = p.Results.MD;
end
if isscalar(p.Results.FA)
        % Input values are min/max number of standard deviations away from
        % the mean
        bounds(5,:) = [-1 1] * p.Results.FA * nanstd(FA(:)) + nanmean(FA(:));
else
    bounds(5,:) = p.Results.FA;
end

% Values that were already NaN in the original dataset.
wasNaN_already = isnan(tensor(:,1,1)); 

% Set values with MD < 1 or MD > 3 or FA < 0.1 to NaN
isOutOfRange = false(n,5);
isOutOfRange(:,1) = L1 < bounds(1,1) | L1 > bounds(1,2);
isOutOfRange(:,2) = L2 < bounds(2,1) | L2 > bounds(2,2);
isOutOfRange(:,3) = L3 < bounds(3,1) | L3 > bounds(3,2);
isOutOfRange(:,4) = MD < bounds(4,1) | MD > bounds(4,2);
isOutOfRange(:,5) = FA < bounds(5,1) | FA > bounds(5,2);

isOutOfAnyRange  = any(isOutOfRange,2);
nOutOfAnyRange   = sum(isOutOfAnyRange);
pctOutOfAnyRange = nOutOfAnyRange/n*100;

% Make values outside the provided range NaNs.
tensor(isOutOfAnyRange,:,:)=NaN;

% Values that are NaN after exclusion.
isNaN_after_excl = isnan(tensor(:,1,1)); 
if nd ==2
    tensor_out = squeeze(tensor);
else
    tensor_out    = reshape(tensor,[sz(1:nd-2) 3 3]);
    isOutOfRange  = reshape(isOutOfRange,[sz(1:nd-2) 5]);
end

% Print results to screen
if p.Results.print == true
    fprintf('in excludeTensor.m:\n')
    fprintf('%d values (%.1f%%) were already NaN \n',...
        sum(wasNaN_already),sum(wasNaN_already) / n*100)
    fprintf('%d values (%.1f%%) were outside of the range\n',...
        nOutOfAnyRange,pctOutOfAnyRange)
    fprintf('%d values (%.1f%%) are NaN after exclusion \n',...
        sum(isNaN_after_excl),sum(isNaN_after_excl) / n*100)
    
end

end

