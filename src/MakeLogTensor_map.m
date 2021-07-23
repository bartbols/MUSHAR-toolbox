function varargout = MakeLogTensor_map( fib_filename,DTI_filename, fa_threshold,varargin )
%MAKELOGTENSOR_MAP Reads the tensor components from the DSI-reconstructed
% .fib-file, applies a Log-Euclidean transform, and saves the log-transformed
% tensor as a nifti-file with the same metadata as the original DTI data.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% May 2021
%
% ----------------- USAGE ----------------- 
% Tensor_filename = MakeLogTensor_map( fib_filename,DTI_filename, fa_threshold)
% or 
% Tensor_filename = MakeLogTensor_map( fib_filename,DTI_filename, fa_threshold,Tensor_filename)
%
% ----------------- INPUT ----------------- 
% - fib_filename     : filename of .fib file as reconstructed with DSI studio
% - DTI_filename     : filename of .nii file with DTI data
% - fa_threshold     : 1x2 vector with lower and upper limit for fractional 
%                      anisotropy values. Values outside of this range are 
%                      set to 0;
%
% Optional 4th input argument:
% - Tensor_filename : filename of new file with eigenvalue data
%
% ----------------- OUTPUT ----------------- 
%-  Tensor_filename   : filename of new file with eigenvalue data
%-  Tensor            : NIfTI structure with eigenvalue data

%% Check inputs
p = inputParser;
addRequired(p,'fib_filename',@(x) contains(x,'.fib'))
addRequired(p,'DTI_filename',@(x) contains(x,'.nii'))
addRequired(p,'fa_threshold',@(x) validateattributes(x,{'numeric'},{'size',[1 2]}))
addOptional(p,'Tensor_filename',[DTI_filename(1:end-7) '_LogTensor.nii.gz'],@(x) contains(x,'.nii'))
parse(p,fib_filename,DTI_filename,fa_threshold,varargin{:});
Tensor_filename = p.Results.Tensor_filename;
%%

% Load the fib file
fprintf('Loading fib-file...')
[path,name,ext] = fileparts(fib_filename);
if strcmp(ext,'.gz')
    gunzip(fib_filename)
    fib_data = load(fib_filename(1:end-3),'-mat');
    delete(fib_filename(1:end-3))
elseif strcmp(ext,'.fib')
    fib_data = load(fib_filename,'-mat');
else
error('Unknown file format %s. The fib-file should have extension .fib or .fib.gz',ext)
end
fprintf(' completed.\n')

fa_map = reshape(fib_data.fa0,fib_data.dimension);

% Load the nifti file with DTI data
DTI_nii = load_untouch_nii(DTI_filename);

% Create new nifti-file with same metadata as original DTI data but with
% the primary eigenvector data (3-channels) as image data
% perm_dim = [2 3 4 1];
% Tensor_map = permute(reshape(fib_data.dir0,[3 fib_data.dimension]),perm_dim);

% Filter the eigenvector map based on FA
txx = reshape(fib_data.txx,fib_data.dimension);
tyy = reshape(fib_data.tyy,fib_data.dimension);
tzz = reshape(fib_data.tzz,fib_data.dimension);
txy = reshape(fib_data.txy,fib_data.dimension);
txz = reshape(fib_data.txz,fib_data.dimension);
tyz = reshape(fib_data.tyz,fib_data.dimension);
txx(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
tyy(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
tzz(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
txy(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
txz(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;
tyz(fa_map < fa_threshold(1) | fa_map > fa_threshold(2)) = 0;

scaling = 1e3;
Tensor_map(:,:,:,1) = txx * scaling;
Tensor_map(:,:,:,2) = tyy * scaling;
Tensor_map(:,:,:,3) = tzz * scaling;
Tensor_map(:,:,:,4) = txy * scaling;
Tensor_map(:,:,:,5) = txz * scaling;
Tensor_map(:,:,:,6) = tyz * scaling;

fprintf('Log transforming the tensor field...')
Tensor_map = logTensor(vec2tensor(Tensor_map));
fprintf(' completed.\n')

% Use all information from the DTI nifti file but overwrite the image data
% with the eigenvalue maps.
% The second dimension needs to be flipped to have a good match in ITK-snap
% between the anatomical data and the eigenvalue map. I'm not sure if this
% also means that the y-direction (or maybe x-direction?) needs to be
% flipped. So it could be that the directions are not stored correctly. For
% segmentation in ITK-snap this doesn't really matter because there will be
% good contrast in colour between tissues with different eigenvectors
% anyway.
% Note 05/02/2018: The eigenvectors are now saved in the ITK coordinate
% system, which is flipped in the x and y coordinates compared to the NIfTI
% coordinate system.

if DTI_nii.hdr.hist.srow_x(1) > 0
    Tensor_map = flip(Tensor_map,1);
end
if DTI_nii.hdr.hist.srow_y(2) > 0
    Tensor_map = flip(Tensor_map,2);
end
    
DTI_nii.img = single(Tensor_map);
DTI_nii.hdr.dime.dim(5) = 6;
DTI_nii.hdr.dime.scl_slope = 1; %/scaling;
DTI_nii.hdr.dime.scl_inter = 0;
DTI_nii.hdr.dime.bitpix = 32;
DTI_nii.hdr.dime.datatype = 16;
save_untouch_nii(DTI_nii,Tensor_filename)

if nargout > 0
    varargout{1} = Tensor_filename;
    if nargout > 1
        varargout{2} = DTI_nii;
    end
end


end % of function

