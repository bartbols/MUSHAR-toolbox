function fa = tensor2fa(tensor)
%TENSOR2FA calculates the fractional anisitropy from a diffusion 
% tensor.
%
% USAGE: fa = tensor2fa(tensor)
%
% INPUT: 
% tensor : n-dimensional tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel.
% 
% OUTPUT:
% fa: (n-2)-dimensional field with fractional anisotropy per voxel.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021

% Get input dimensions.
nd = ndims(tensor);
sz = size(tensor);

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-2));
tensor = reshape(tensor,[n,3,3]);

fa = NaN(n,1);
for i = 1 : n
    
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        continue
    end
    % Decompose into rotation matrix v and eigenvalue matrix l.
    [~,D] = eig(S);
    L = diag(D);
    if any(L <= 0)
        continue
    end
    fa(i) = sqrt(1/2) * sqrt((L(1) - L(2)).^2 + (L(2) - L(3)).^2 + (L(3) - L(1)).^2 ) ./ ...
                         sqrt(L(1).^2 + L(2).^2 + L(3).^2);
end

% Put tensor field back in original (input) dimensions.
if nd == 3
    fa = squeeze(fa);
else
    fa    = reshape(fa,sz(1:nd-2));
end


end

