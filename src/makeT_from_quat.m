function [ T,R ] = makeT_from_quat( img )
%MAKET_FROM_QUAT builds the 4x4 spatial transformation matrix from the
%header information in the NIfTI structure 'img'.

if ~isstruct(img)
    img = load_untouch_nii(img);
end

b = img.hdr.hist.quatern_b;
c = img.hdr.hist.quatern_c;
d = img.hdr.hist.quatern_d;
a = sqrt(1.0-(b*b+c*c+d*d));
a = real(a);

R =[ a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c;...
     2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b;...
     2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b];
 
% Scale R by voxel dimensions and put in T
T = eye(4);
T(1:3,1:3) = R * diag(img.hdr.dime.pixdim(2:4));

% flip the z-direction if qfac = 1;
qfac = img.hdr.dime.pixdim(1);
T(:,3) = qfac * T(:,3);

% Add the offsets
T(1:3,4) = [img.hdr.hist.qoffset_x;...
            img.hdr.hist.qoffset_y;...
            img.hdr.hist.qoffset_z];

end

