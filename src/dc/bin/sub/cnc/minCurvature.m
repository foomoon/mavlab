function r = minCurvature(X,Y,Z);
[g,m] = surfature(X,Y,Z);
m = m(:);
m(isnan(m)) = [];
m = max(m);
r = 1/m;