function [x,y,z,v,s] = offsetSurface(X,Y,Z,r)
% OFFSETSURFACE - Creates an offset surface
%   [x,y,z] = offsetSurface(X,Y,Z,R) will generate a new surface from the
%   surface X,Y,Z.  The new surface is offset in the direction normal to
%   the surface by distance R. 

global dcDebug

[K,H] = surfature(X,Y,Z);
R = abs(1./H);
Violations = R < r;
ii = find(Violations);

% [nx,ny,nz] = surfnorm(X,Y,Z);
[nx,ny,nz] = normals(X,Y,Z);
x = X - nx*r;
y = Y - ny*r;
z = Z - nz*r;

f = zeros(size(nx));
for i=1:length(ii),
    f = f + (sqrt((X-X(ii(i))).^2 + (Y-Y(ii(i))).^2) < r) ;
end
allViolations = logical(f);
s.x = X; s.x(~allViolations) = nan;
s.y = Y; s.y(~allViolations) = nan;
s.z = Z; s.z(~allViolations) = nan;
h = H; h(~allViolations) = nan;
h = unique(sort(round(10000*abs(1./h(~isnan(h))))/10000));
s.r = h(h<r);
v = allViolations;

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  Warning: %0.0f surface violation(s)\n',sum(v(:)))
end