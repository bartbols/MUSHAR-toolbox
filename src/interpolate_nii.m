function [ values ] = interpolate_nii( nii,points,method )
%INTERPOLATE_NII interpolate the NIfTI image 'nii' at the n x 3 coordinates
%provided by 'points'). Points should be provided in global coordinates
%(not voxel coordinates)
%
% values = interpolate_nii(nii,points)
% values = interpolate_nii(nii,points,method)
%
% If no third input argument is provided, method = 'linear'.
%
% Bart Bolsterlee
% February 2018
% Neuroscience Research Australia
%
% Change log:
% 18/10/2018, BB: Added interpolation method as third, optional input.

if nargin < 3
    method = 'linear';
    
end
[m,n] = size(points);
if m == 3 && n ~=3
    % Transpose
    points = points';
    TransposeFlag = true;
else
    TransposeFlag = false;
end
    
if ~isstruct(nii)
    % If nii is not a structure, it is assumed to be a filename of a NIfTI
    % file. The NIfTI file is now loaded.
    nii = load_untouch_nii(nii);
end
% Get the transformation matrix from the NIfTI header.
% if all(nii.hdr.hist.srow_x == 0)
T = makeT_from_quat(nii);
% else
%     T = [nii.hdr.hist.srow_x;...
%          nii.hdr.hist.srow_y;...
%          nii.hdr.hist.srow_z;...
%          0 0 0 1];
% end

% Convert points to voxel coordinates
voxel_coords = [points ones(size(points,1),1)] / (T');

% Interpolate the image at these voxel locations
I = squeeze(double(nii.img));
imdim = size(I);

if numel(imdim) == 3
    values = interp3(I,voxel_coords(:,2)+1,voxel_coords(:,1)+1,voxel_coords(:,3)+1,method);
elseif numel(imdim) == 4
    values = zeros(size(points,1),imdim(4));
    for i = 1 : imdim(4)
        values(:,i) = interp3(I(:,:,:,i),voxel_coords(:,2)+1,voxel_coords(:,1)+1,voxel_coords(:,3)+1,method);
    end
else
    error('Only 3D and 4D images are supported.')
end

if TransposeFlag == true
    values = values';
end

end

