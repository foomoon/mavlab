function ff = formfactor(tr,formtype),
% FORMFACTOR Computation Skin Friction Coeficient
%   FORMFACTOR returns the value at which drag needs to be scaled by to 
%   account for the frontal profile of the object.
%
%   FORMFACTOR(TR,formtype) returns the form factor scaling coefficient,
%   where TR is the thickness ratio (t/c or d/c) of the object and formtype
%   is the type of object:
%
%       formtype
%       --------
%       'w' - Wing-like objects
%       'f' - Fuselage-like objects   
%
%   Input:
%       TR:           (scaler) Thickness Ratio
%       Formtype:     (string) Formfactor Type
%
%   Output:
%       FF:           (scaler) Form Factor Scaling Coefficient
%
%   See also FRICTION

if strcmp(lower(formtype(1)),'w'),
    ff = 1.0 + 1.8*tr + 50*tr^4;
elseif strcmp(lower(formtype(1)),'f'),
    ff = 1.0 + 1.5*tr^1.5 + 50*tr^3;
else
    ff = 0;
end