function [Xdeg] = degrees(Xrad)
% DEGREES Convert from radians to degrees
%   [Xdeg] = degrees(Xrad)
%   converts a value in radians (Xrad) to a value in degrees (Xdeg)
%
%   see also RADIANS, VISUALANGLE

Xdeg = Xrad * (180/pi);