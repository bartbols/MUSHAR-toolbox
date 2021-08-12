function vec = tensor2vec(tensor)
%TENSOR2VEC Reformats the tensor field in 3x3 full tensor format to tensor 
% field in 1x6 vector format.
% 
% USAGE:
% tensor = vec2tensor(vec)
% INPUT:
% tensor: (i,j,k) x 3x3 array with tensor field in full tensor format.
% OUTPUT:
% vec: (i,j,k) x 6 array with tensor field in vector format: [txx tyy tzz txy txz tyz]

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
tensor = reshape(tensor,[n,9]);

vec = tensor(:,[1 5 9 2 3 6]);
if nd > 2
    vec = reshape(vec,[sz(1:nd-2) 6]);
end

end

