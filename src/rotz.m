  function [Rz]=rotz(th)
% Program that creates a rotation matrix for a rotation round the z Axis
Rz(1,1)=cos(th);
Rz(1,2)=-sin(th);
Rz(2,1)= sin(th);
Rz(2,2)= cos(th);
Rz(3,3)=1;
