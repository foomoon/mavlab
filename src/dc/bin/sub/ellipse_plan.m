function [x,y] = ellipse_plan(a,b),
% ELLIPSE_PLAN Create an Ellipse
%   ELLIPSE_PLAN Calculate curvature of an ellipse given the semi-major and
%   semi-minor axes.
%
%   [x,y] = ELLIPSE_PLAN(a,b)
%
%   Input
%       a:   (scaler) semi-major axis
%       b:   (scaler) semi-minor axis
%
%   Output
%       x:   (vector) x-coordinates
%       y:   (vector) y-coordinates
%
% See also: MAVLAB

pts=41;
theta = 0:(pi/(pts-1)):pi;
x = -a*(cos(theta));
y = (b)*sqrt((1-(x.^2)/(a^2)));

