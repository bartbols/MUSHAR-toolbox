function surface_out = transformSurface(surface_in,T)
%TRANSFORMSURFACE Transforms the vertices of the input surface with the
%provided transformation matrix.
%
% USAGE:
% surface_out = transformSurface(surface_in,T)
% 
% INPUT:
% surface_in : stl-filename OR a structure with fields 'Points' and 'ConnectivityList'
% T          : 4x4 transformation matrix.
%
% OUTPUT:
% surface_out : structure with fields 'Points' and 'ConnectivityList'
%               with the transformed surface.
%               The input points are transformed following:
%               points_out = (points_in - T(1:3,4)') * T(1:3,1:3)
%
% Bart Bolsterlee, Neuroscience Research Australia
% July 2021
%              

if ~isstruct(surface_in)
    surface_in = stlread(surface_in);
end

surface_out.Points = (surface_in.Points - T(1:3,4)') * T(1:3,1:3);
surface_out.ConnectivityList = surface_in.ConnectivityList;

end

