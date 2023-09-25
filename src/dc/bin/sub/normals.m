function [nx,ny,nz] = normals(x,y,z)

[nx,ny,nz] = quadcross(x,y,z);
[nx,ny,nz] = unitvec(nx,ny,nz);
nx = padmat(nx);
ny = padmat(ny);
nz = padmat(nz);
nx = quadavg(nx);
ny = quadavg(ny);
nz = quadavg(nz);

function xavg = quadavg(x)
[m,n] = size(x);
i = 1:m-1;
j = 1:n-1;

xavg = (x(i,j) + x(i+1,j) + x(i,j+1) + x(i+1,j+1) ) / 4;


function [nx,ny,nz] = quadcross(x,y,z)
[m,n] = size(x);
i = 1:m-1;
j = 1:n-1;
ax = (x(i,j)-x(i+1,j+1));
bx = (x(i+1,j)-x(i,j+1));

ay = (y(i,j)-y(i+1,j+1));
by = (y(i+1,j)-y(i,j+1));

az = (z(i,j)-z(i+1,j+1));
bz = (z(i+1,j)-z(i,j+1));

nx = -(ay.*bz - az.*by);
ny =  (ax.*bz - az.*bx);
nz = -(ax.*by - ay.*bx);

function [x,y,z] = unitvec(x,y,z)
mag = sqrt(x.^2 + y.^2 + z.^2);
mag = max(mag,eps);
x = x./mag;
y = y./mag;
z = z./mag;

function x = padmat(x)
x = [x(:,1) x x(:,end)];
x = [x(1,:); x; x(end,:)];