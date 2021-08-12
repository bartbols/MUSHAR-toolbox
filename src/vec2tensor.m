function tensor = vec2tensor(vec)
%VEC2TENSOR Reformats the tensor field in 1x6 vector format to a full 3x3 
% tensor. It expects the vector to be in the format: [txx tyy tzz txy txz tyz]
%
% USAGE:
% tensor = vec2tensor(vec)
% INPUT:
% vec: (i,j,k) x 6 array with tensor field in vector format
% OUTPUT:
% tensor: (i,j,k) x 3x3 array with tensor field in full tensor format.
%
% Bart Bolsterlee
% Neuroscience Research Australia
% May 2021

% Get input dimensions.
nd = ndims(vec);
sz = size(vec);

if sz(nd) ~= 6
    error('Input array should have 6 as last dimension')
end

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-1));
vec = reshape(vec,[n,6]);

tensor = vec(:,[1 4 5 4 2 6 5 6 3]);
tensor = reshape(tensor,[sz(1:nd-1) 3 3]);
if nd==2
    tensor = squeeze(tensor);
end

end

