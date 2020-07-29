function [Xrad] = radians (Xdeg)
% RADIANS Convert from degrees to radians
%  [Xrad] = radians (Xdeg)
%   converts a value in degrees (Xdeg) to a value in radians (Xrad)
%
%   see also DEGREES, VISUALANGLE

Xrad = Xdeg * (pi/180);