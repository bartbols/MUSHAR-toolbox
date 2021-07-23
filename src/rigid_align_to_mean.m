function [Xt,tform] = rigid_align_to_mean(X)
%RIGID_ALIGN_TO_MEAN rigidly aligns all shapes in array X to the mean
%shape. X is assumed to be an array of point-to-point corresponding points
%on a surface.

% % Align all shapes to the mean using Procrustes analysis
X_mean = reshape(mean(X,2),[],3);
s=size(X,2);
Xt = NaN(size(X));
for nr = 1 : s
    [d,X_aligned,tform(nr)] = procrustes(...
        X_mean,...
        reshape(X(:,nr),[],3),...
        'scaling',false,...
        'reflection',false);
    Xt(:,nr) = X_aligned(:);
end

end % of function

