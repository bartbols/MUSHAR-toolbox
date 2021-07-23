function [ev1,ev2,ev3] = tensor2ev(tensor)
%TENSOR2EV calculates the eigenvectors from a diffusion tensor field.
%
% USAGE: [ev1,ev2,ev3] = tensor2ev(tensor)
%
% INPUT: 
% tensor : n-dimensional tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel.
% 
% OUTPUT:
% ev1 : (n-1)-dimensional field with in the last dimension the 3 components
%       of the primary eigenvector.
% ev2 : (n-1)-dimensional field with in the last dimension the 3 components
%       of the secondary eigenvector.
% ev3 : (n-1)-dimensional field with in the last dimension the 3 components
%       of the tertiary eigenvector.
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

ev1 = NaN(n,3);
ev2 = NaN(n,3);
ev3 = NaN(n,3);
for i = 1 : n
%     if mod(i,100)==0
%         fprintf('%d of %d\n',i,n)
%     end
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        continue
    end
    
    % Decompose into rotation matrix v and eigenvalue matrix l.
    [v,l] = eig(S);
    
    if any(diag(l) <= 0)
        % If one of the eigenvalues is negative, the eigenvectors will be set
        % to NaN for that voxel.
        continue
    end    
    
    % Get primary eigenvector.
    [~,idx] = sort(diag(l),'descend');
    ev1(i,:) = v(:,idx(1))';
    ev2(i,:) = v(:,idx(2))';
    ev3(i,:) = v(:,idx(3))';
end

% Put tensor field back in original (input) dimensions.
if nd > 2
    ev1    = reshape(ev1,[sz(1:nd-2) 3]);
    ev2    = reshape(ev2,[sz(1:nd-2) 3]);
    ev3    = reshape(ev3,[sz(1:nd-2) 3]);
end

end