function log_vec = logTensor(tensor)
%LOGTENSOR calculates the Log-Euclidean transformation of a diffusion 
% tensor following procedures described in: 
%
% Arsigny, V., Fillard, P., Pennec, X., Ayache, N., 2006. Log-Euclidean
% metrics for fast and simple calculus on diffusion tensors. Magn Reson Med
% 56 (2),p 411-421.http://doi.org/10.1002/mrm.20965
%
% USAGE: log_vec = logTensor(tensor)
%
% INPUT: 
% tensor : n-dimensional tensor field. The last two dimensions should be 3x3,
%     containing the tensor for that voxel.
% 
% OUTPUT:
% log_vec: (n-1) array with the last dimension the 6 elements of 
%          Log-Euclidean tensor in vector representation.
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

log_tensor = NaN(size(tensor));
log_vec    = NaN(n,6);
for i = 1 : n
    
    % Step 1: decompose into rotation matrix R and eigenvalue matrix D
    % (eigenvalues on diagonal).
    S = squeeze(tensor(i,:,:));
    if any(isnan(S(:)))
        % Set value to NaN if any of the eigenvalues if smaller than zero
        % or
        log_tensor(i,:,:) = NaN;
        log_vec(i,1:6) = NaN;
        continue
    end    
    
    [R,D] = eig(S);
    
    if any(diag(D) <= 0)
        % Set value to NaN if any of the eigenvalues if smaller than zero
        % or
        log_tensor(i,:,:) = NaN;
        log_vec(i,1:6) = NaN;
        continue
    end
    % Step 2: log-transform eigenvalues
    Dstar = diag(log(diag(D)));
    
    % Step 3: recompose tensor
    log_tensor(i,:,:) = R*Dstar*R'; % note that R in Arsigny et al. equals R' here
    %                         % because of Matlab's convention of representing
    %                         % the eigenvector matrix in column format.
    
    % Step 4: put in vector form
    log_vec(i,1:6) = [...
        log_tensor(i,1,1)...
        log_tensor(i,2,2)...
        log_tensor(i,3,3),...
        sqrt(2)*log_tensor(i,1,2),...
        sqrt(2)*log_tensor(i,1,3),...
        sqrt(2)*log_tensor(i,2,3)];
end

% Put tensor field back in original (input) dimensions.
log_vec    = reshape(log_vec,[sz(1:nd-2),6]);

end

