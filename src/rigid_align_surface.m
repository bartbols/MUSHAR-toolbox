function [FV_aligned,T] = rigid_align_surface(refSurface,targetSurface,alignedSurface)
%RIGID_ALIGN_SURFACE dalls ShapeWorks to rigidly align a reference to a
%target surface using the iterative closest point algorithm. The aligned
%surface and the rotation matrix is recovered.
%
% USAGE
% [FV_aligned,T] = rigid_align_surface(refSurface,targetSurface,alignedModel)
% 
% INPUT:
% refSurface     : stl filename of the reference surface
% targetSurface  : stl filename of the target surface
% alignedSurface : stl filename of the aligned surface (will be created)
% 
% OUTPUT         
% FV_aligned     : structure with Points and ConnectivityList of aligned
%                  surface
% T              : 4x4 transformation matrix that rigidly aligns the
%                  reference with the target surface


% alignedSurface = fullfile(tempdir,'aligned.stl');
if exist(alignedSurface,'file')==2;delete(alignedSurface);end
shapeworks_cmd = sprintf('shapeworks read-mesh --name=%s transform-mesh --target=%s --type=%s --method=icp --iterations=%d write-mesh --name=%s',...
        targetSurface,refSurface,'rigid',100,alignedSurface);
[status,cmdout] = system(shapeworks_cmd);

if contains(cmdout,'is not recognized as an internal or external command')
    error('shapeworks not found as an external commant. Please install ShapeWorks and try again.')
end

% Recover the transformation matrix 
FV_aligned = stlread(alignedSurface);
FV_target  = stlread(targetSurface);
[T,eps] = estimateRigidTransform(FV_target.Points',FV_aligned.Points');
% stlwrite(FV_aligned,alignedModel);

