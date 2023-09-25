function [ac,mean_ac] = mac(xl,yl,xt,yt),
% MAC Mean Aerodynamc Center
%   MAC calculates the aerodynamic center of a wing given the leading and
%   trailing edge coordinates.
%
%   [ac,mean_ac] = MAC(xl,yl,xt,yt)
%
%   Input:
%       xl: (vector) Leading edge x coordinates
%       yl: (vector) Leading edge y coordinates
%       xt: (vector) Trailing edge x coordinates
%       yt: (vector) Trailing edge y coordinates
%
%   Output:
%       ac: (scaler) Aerodynamic Center from (0,0,0)
%
% See also: MAVLAB

% c = yl-yt;
% cx = yl.*c;
% c2 = c.^2;
% 
% xa = sum(c2)/sum(c);
% xb = sum(cx)/sum(c);
% 
% out = xb - (xa/4);

c = yl-yt;
cx = yl.*c;
c2 = c.^2;

mean_ac = sum(c2)/sum(c);
xb = sum(cx)/sum(c);
ac = xb - (mean_ac/4);
