  function [Ry]=roty(th)
% Program that creates a rotation matrix for a rotation round the y Axis
Ry(1,1)=cos(th);
Ry(1,3)=sin(th);
Ry(2,2)=1;
Ry(3,1)=-sin(th);
Ry(3,3)= cos(th);
