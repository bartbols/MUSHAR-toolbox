function [l1,l2,l3,md,fa] = tensor2diffprop(tensor)
%TENSOR2DIFFPROP calculates the diffusion properties from a diffusion tensor
%field. Diffusion properties are the primary, secondary and tertiary
%eigenvalues, the mean diffusivity and the fractional anisotropy.
%
% USAGE: [l1,l2,l3,md,fa] = tensor2diffprop(tensor)
%
% INPUT: 
% tensor : n-dimensional tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel.
% 
% OUTPUT:
% lambda1 : (n-2)-dimensional field with the primary eigenvalue per voxel.
% lambda2 : (n-2)-dimensional field with the secondary eigenvalue per voxel.
% lambda3 : (n-2)-dimensional field with the tertiary eigenvalue per voxel.
% lambda3 : (n-2)-dimensional field with the mean diffusivity per voxel.
% fa      : (n-2)-dimensional field with the fractional anisotropy per voxel.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021

% Get input dimensions.
nd = ndims(tensor);
sz = size(tensor);

if any(sz((nd-1:nd)) ~= [3 3])
    error('Last two dimensions should be 3x3')
end

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-2));
tensor = reshape(tensor,[n,3,3]);

l1 = NaN(n,1);
l2 = NaN(n,1);
l3 = NaN(n,1);
for i = 1 : n
%     if mod(i,100)==0
%         fprintf('%d of %d\n',i,n)
%     end
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        continue
    end
    
    % Decompose into rotation matrix v and eigenvalue matrix l.
    [~,l] = eig(S);
    
    if any(diag(l) <= 0)
        % If one of the eigenvalues is negative, the eigenvalues will be set
        % to NaN for that voxel.
        continue
    end    
    
    % Get primary eigenvector.
    l = sort(diag(l),'descend');
    l1(i) = l(1);
    l2(i) = l(2);
    l3(i) = l(3);
end

                
% Put eigenvalue fields back in original (input) dimensions.
if nd > 3
    l1    = reshape(l1,sz(1:nd-2));
    l2    = reshape(l2,sz(1:nd-2));
    l3    = reshape(l3,sz(1:nd-2));
end

md = (l1 + l2 + l3)/3;
fa = sqrt(1/2) * sqrt((l1 - l2).^2 + (l2 - l3).^2 + (l3 - l1).^2 ) ./ sqrt(l1.^2 + l2.^2 + l3.^2);


end

