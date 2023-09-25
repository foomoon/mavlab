function y = dihedral(x,theta,polyhedral)
% DIHEDRAL Calculate dihedral angle of wing
%     DIHEDRAL(x,theta,poly) returns the spanwise elevation of a wing with
%     respect to an input dihedral angle (theta) and polyhedral ratio
%     (poly).
%
%     Z = DIHEDRAL(x,theta,poly)
%
%     Input:
%         x:       (vector) spanwise coordinates
%         theta:   (scaler) value of dihedral angle (degrees)
%         poly:    (scaler) ratio of the span at which dihedral starts
%                (0 < poly < 1)
%
%     Output:
%         Z:    (vector) dihedral values
%
% See also: MAVLAB

theta=theta*pi/180;
for i=1:length(x)
    y(i)= tan(theta)*(x(i)-polyhedral) ; %( abs(x(i)) - polyhedral*abs(x(1)) )*tan(theta);
    if abs(x(i))<polyhedral
        y(i)=0;
    end
end

