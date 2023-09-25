function x = tol2res(tol,r),

% r = D/2;

% r-tol = r*cos(asin(x/(2*r)))
% 2*r*sin(acos(1-tol/r)) = (x)

x = 2*r*sin(acos(1-tol/r));