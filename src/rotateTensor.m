function tensor_rot = rotateTensor(tensor,R)
%ROTATETENSOR Rotates tensor field with rotation matrix R following:
% Trot = R' * tensor * R
%
% INPUT:
% tensor field

%% Get input dimensions.
nd = ndims(tensor);
sz = size(tensor);

isvec=false;
if sz(nd) == 6
    % Input is in vector format. Convert to full tensor format for
    % processing.
    isvec=true;
    tensor = vec2tensor(tensor);
    nd=nd+1;
end

% Reshape tensor field into one long array of tensors for easy processing.
n = prod(sz(1:nd-2));
tensor = reshape(tensor,[n,3,3]);

tensor_rot = NaN(size(tensor));
tensor_rot(:,1,1) = R(1,1)*(R(1,1)*tensor(:,1,1) + R(2,1)*tensor(:,1,2) + R(3,1)*tensor(:,1,3)) + R(2,1)*(R(1,1)*tensor(:,2,1) + R(2,1)*tensor(:,2,2) + R(3,1)*tensor(:,2,3)) + R(3,1)*(R(1,1)*tensor(:,3,1) + R(2,1)*tensor(:,3,2) + R(3,1)*tensor(:,3,3));
tensor_rot(:,1,2) = R(1,1)*(R(1,2)*tensor(:,1,1) + R(2,2)*tensor(:,1,2) + R(3,2)*tensor(:,1,3)) + R(2,1)*(R(1,2)*tensor(:,2,1) + R(2,2)*tensor(:,2,2) + R(3,2)*tensor(:,2,3)) + R(3,1)*(R(1,2)*tensor(:,3,1) + R(2,2)*tensor(:,3,2) + R(3,2)*tensor(:,3,3));
tensor_rot(:,1,3) = R(1,1)*(R(1,3)*tensor(:,1,1) + R(2,3)*tensor(:,1,2) + R(3,3)*tensor(:,1,3)) + R(2,1)*(R(1,3)*tensor(:,2,1) + R(2,3)*tensor(:,2,2) + R(3,3)*tensor(:,2,3)) + R(3,1)*(R(1,3)*tensor(:,3,1) + R(2,3)*tensor(:,3,2) + R(3,3)*tensor(:,3,3));
tensor_rot(:,2,1) = R(1,2)*(R(1,1)*tensor(:,1,1) + R(2,1)*tensor(:,1,2) + R(3,1)*tensor(:,1,3)) + R(2,2)*(R(1,1)*tensor(:,2,1) + R(2,1)*tensor(:,2,2) + R(3,1)*tensor(:,2,3)) + R(3,2)*(R(1,1)*tensor(:,3,1) + R(2,1)*tensor(:,3,2) + R(3,1)*tensor(:,3,3));
tensor_rot(:,2,2) = R(1,2)*(R(1,2)*tensor(:,1,1) + R(2,2)*tensor(:,1,2) + R(3,2)*tensor(:,1,3)) + R(2,2)*(R(1,2)*tensor(:,2,1) + R(2,2)*tensor(:,2,2) + R(3,2)*tensor(:,2,3)) + R(3,2)*(R(1,2)*tensor(:,3,1) + R(2,2)*tensor(:,3,2) + R(3,2)*tensor(:,3,3));
tensor_rot(:,2,3) = R(1,2)*(R(1,3)*tensor(:,1,1) + R(2,3)*tensor(:,1,2) + R(3,3)*tensor(:,1,3)) + R(2,2)*(R(1,3)*tensor(:,2,1) + R(2,3)*tensor(:,2,2) + R(3,3)*tensor(:,2,3)) + R(3,2)*(R(1,3)*tensor(:,3,1) + R(2,3)*tensor(:,3,2) + R(3,3)*tensor(:,3,3));
tensor_rot(:,3,1) = R(1,3)*(R(1,1)*tensor(:,1,1) + R(2,1)*tensor(:,1,2) + R(3,1)*tensor(:,1,3)) + R(2,3)*(R(1,1)*tensor(:,2,1) + R(2,1)*tensor(:,2,2) + R(3,1)*tensor(:,2,3)) + R(3,3)*(R(1,1)*tensor(:,3,1) + R(2,1)*tensor(:,3,2) + R(3,1)*tensor(:,3,3));
tensor_rot(:,3,2) = R(1,3)*(R(1,2)*tensor(:,1,1) + R(2,2)*tensor(:,1,2) + R(3,2)*tensor(:,1,3)) + R(2,3)*(R(1,2)*tensor(:,2,1) + R(2,2)*tensor(:,2,2) + R(3,2)*tensor(:,2,3)) + R(3,3)*(R(1,2)*tensor(:,3,1) + R(2,2)*tensor(:,3,2) + R(3,2)*tensor(:,3,3));
tensor_rot(:,3,3) = R(1,3)*(R(1,3)*tensor(:,1,1) + R(2,3)*tensor(:,1,2) + R(3,3)*tensor(:,1,3)) + R(2,3)*(R(1,3)*tensor(:,2,1) + R(2,3)*tensor(:,2,2) + R(3,3)*tensor(:,2,3)) + R(3,3)*(R(1,3)*tensor(:,3,1) + R(2,3)*tensor(:,3,2) + R(3,3)*tensor(:,3,3));

% % Code to check that the vectorized version is correct.
% tensor_rot2 = NaN(size(tensor));
% for i =1 : n
%     tensor_rot2(i,1:3,1:3) = R'*squeeze(tensor(i,1:3,1:3))*R;
% end
% isequal(tensor_rot2,tensor_rot) % Should be equal (within a very narrow
%                                 % margin).
% max(abs(tensor_rot(:)-tensor_rot2(:)))
% 
% Put tensor field back in original (input) dimensions.
tensor_rot = reshape(tensor_rot,[sz(1:nd-2),3,3]);
if isvec == true
    % Convert back to vector format.
    tensor_rot = tensor2vec(tensor_rot);
end

end

