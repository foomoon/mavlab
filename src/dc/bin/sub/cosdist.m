function x = cosdist(num_pts),
% COSDIST Cosine distribution of points
%   COSDIST(num_pts) Builds an array of points between 0-1 with a cosine
%   distribution.  This function is similar to linspace, but does not 
%   require upper and lower limits.
%
%   X = COSDIST(41) produces a (1 x 41) vector of points between 0 and 1.
%
%   Input:
%       num_pts:  (scaler) number greater than 0
%
%   Output: 
%       X:  (vector) [1 x num_pts] distribution of points
%
% See also: MAVLAB

angle_incriment = (pi)/(num_pts-1);
i=1;

for theta=0:angle_incriment:pi,
    x(i) = -.5*(cos(theta)-1);
    i = i + 1;
end