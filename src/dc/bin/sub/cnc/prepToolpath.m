function [X,Y,Z] = prepToolpath(X,Y,Z,step,res)

global dcDebug

% Find Repeated points
% 1st derivative = 0
dx = diff(X(:));
dy = diff(Y(:));
dz = diff(Z(:));

% In stead of using the expression dx == 0, etc.  It is more appropriate
% to asses weather the difference is less than the machine precision.  If
% so, then the derivative equivalently zero. The term "eps" is not a
% variable, it is a function provided by MATLAB. Look at the help for eps
% for more information.
rep = (abs(dx) < eps & abs(dy) < eps & abs(dz) < eps);
rep = [false; rep];

% Find contiguous points lying on the same line
% 2nd derivative = 0
ddx = diff(dx);
ddy = diff(dy);
ddz = diff(dz);

% Same as above, we use abs(ddx) < eps in stead of ddx == 0
con = (abs(ddx) < eps & abs(ddy) < eps & abs(ddz) < eps);
con = [false; con; false];

% Find nan's
nn = isnan(X) | isnan(Y) | isnan(Z);

% Remove Selected Points
rem = rep | con | nn;
X(rem) = [];
Y(rem) = [];
Z(rem) = [];

% Find any places where the distance between points is much larger than the res
D = sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2);
n = 10;
rem = D > n*step & D > n*res;
if rem(end)
    e = true;
else
    e = false;
end
gap = [rem; e];
X(gap) = [];
Y(gap) = [];
Z(gap) = [];

% Round off significant figures
n = 4;
X = round(X*10^n)/(10^n);
Y = round(Y*10^n)/(10^n);
Z = round(Z*10^n)/(10^n);

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  Removed %0.0f repeated point(s)\n',sum(rep))
    fprintf(1,'  Removed %0.0f linear point(s)\n',sum(con))
    fprintf(1,'  Removed %0.0f NaN point(s)\n',sum(nn))
    fprintf(1,'  Removed %0.0f poorly spaced point(s)\n',sum(gap))
    fprintf(1,'  Rounded to %0.0f decimal places\n',n)
end



