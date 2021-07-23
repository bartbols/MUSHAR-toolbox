  function [Rx]=rotx(th)
% Program that creates a rotation matrix for a rotation round the x Axis
Rx(2,2)=cos(th);
Rx(2,3)=-sin(th);
Rx(3,2)= sin(th);
Rx(3,3)= cos(th);
Rx(1,1)=1;
