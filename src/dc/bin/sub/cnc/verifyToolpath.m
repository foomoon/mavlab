function [h] = verifyToolpath(X,Y,Z,D,varargin),
% VERIFYTOOLPATH - Simulation of Toolpath
%   VERIFYTOOLPATH is a Computer Numerically Controlled (CNC) toolpath
%   verification simulation.  It attempts to render the actual milling
%   process of a surface using a hemispherical cutting tool known as a
%   ball-endmill.
%   VERIFYTOOLPATH(X,Y,Z,D) requires the toolpath given by vectors X,Y,Z,
%   and the diameter of the ball-endmill (D).
%   VERIFYTOOLPATH(...,REFRESH), specifies how often the rendering will be
%   updated during milling.  This parameter's default value is 0.01, and
%   ranges from 0-1, where 0 will update continously, and 1 will update
%   only when calculations are finished.
%
%


% Create rectangular Workpiece based on maximum dimensions and tool
% diameter (D)
minX = min(X(:)) - D/2;
maxX = max(X(:)) + D/2;
minY = min(Y(:)) - D/2;
maxY = max(Y(:)) + D/2;
minZ = min(Z(:)) - D/2;

% Resolution of model 
ptsPerInch   = 75;

xGridDensity = round((maxX-minX)*ptsPerInch);
yGridDensity = round((maxY-minY)*ptsPerInch);

xGrid = linspace(minX,maxX,xGridDensity);
yGrid = linspace(minY,maxY,yGridDensity);

[xBlock,yBlock] = meshgrid(xGrid,yGrid);
[zBlock] = zeros(size(xBlock));

createBottom(minX,maxX,minY,maxY,minZ)

% Initialize Workpiece model
h = surf(xBlock,yBlock,zBlock,...
    'edgecolor','none',...
    'facecolor',[.2 .8 1],...
    'facelight','gouraud',...
    'specularstrength',0,...
    'ambientstrength',.3);
light('position',[.3 .3 .3]);
axis equal
view(38,42);
set(gca,'drawmode','fast');

drawnow

%% Create ball endmill
% [xb,yb,zb,m,n] = endMill(D);

N = length(X);

%% Plotting Settings
update = 1;
if nargin == 4,    
    upFreq = round(.01*N);
else
    upFreq = round(N*varargin{1});
    if isempty(upFreq),
        update = 0;
    end
end

R = D/2;
Path = [X(:) Y(:) Z(:)];

%% Solve
maxiter = length(X) - 1;
for i=1:maxiter,
    pc = i:i+1;
    [zTool,yRef,xRef] = traverse(R,xBlock,yBlock,Path(pc,:));
    zBlock(yRef,xRef) = min(zBlock(yRef,xRef),zTool);
    if update & ~mod(i,upFreq),
        set(h,'Zdata',(zBlock));
        drawnow
    end
    workbar(i/maxiter,'Computing Toolpath');
end
set(h,'Zdata',zBlock);
drawnow





function [xb,yb,zb,mm,nn] = endMill(D),
[xb,yb,zb] = sphere(20); 
mm = round(size(xb,1)/4); 
nn = round(size(xb,2)/2); 

m = round(size(xb,1)/2); 
n = round(size(xb,2)/1);

xb = xb(1:m,1:n)*D/2;
yb = yb(1:m,1:n)*D/2;
zb = zb(1:m,1:n)*D/2 + D/2;



function [Z,yRef,xRef] = traverse(R,xBlock,yBlock,pos),

% Point i and point i+1
x = pos(:,1);
y = pos(:,2);
z = pos(:,3);

% Find unit vector from point i to point i+1
u = diff(pos);
mag = sqrt(u(1)^2 + u(2)^2 + u(3)^2);
mag = max(eps,mag);
u = u/mag;

% Create Extents of cutting region
maxX = x + R;  maxX = max(maxX);
minX = x - R;  minX = min(minX);
maxY = y + R;  maxY = max(maxY);
minY = y - R;  minY = min(minY);

% Find local xBlock & yBlock region enclosing local cutting area
xRef = xBlock(1,:) >= minX & xBlock(1,:) <= maxX;
yRef = yBlock(:,1) >= minY & yBlock(:,1) <= maxY;
% Local cutting region
X = xBlock(yRef,xRef);
Y = yBlock(yRef,xRef);

% Compute semi-sphere at local grid (P1)
Z1 = -sqrt(-(X-x(1)).^2 - (Y-y(1)).^2 + R^2) + z(1);
% Compute semi-sphere at local grid (P2)
Z2 = -sqrt(-(X-x(2)).^2 - (Y-y(2)).^2 + R^2) + z(2);
% Z component of Cylindar about arbitray axis u
Zc =  (-u(2).^2.*X.*u(1).*u(3)+z(1)-X.*u(1).^3.*u(3)+u(1).^3.*x(1).*u(3)+u(3).^3.*u(2).*y(1)-2.*u(3).*u(1).*x(1)-2.*u(3).*u(2).*y(1)-u(1).^2.*Y.*u(2).*u(3)+u(1).^2.*u(3).^2.*z(1)+2.*u(3).*X.*u(1)-Y.*u(2).^3.*u(3)+u(2).^3.*y(1).*u(3)+u(2).^2.*u(3).^2.*z(1)+2.*u(3).*Y.*u(2)-u(3).^3.*X.*u(1)+u(3).^3.*u(1).*x(1)-u(3).^3.*Y.*u(2)+u(1).^2.*u(2).*y(1).*u(3)+u(2).^2.*u(1).*x(1).*u(3)+u(3).^4.*z(1)-2.*u(3).^2.*z(1)-(-4.*u(3).^2.*Y.*y(1)-u(2).^2.*X.^2.*u(1).^2-u(1).^2.*u(2).^2.*y(1).^2-u(2).^2.*u(3).^2.*X.^2-4.*u(3).^2.*X.*x(1)-u(3).^4.*Y.^2+2.*u(3).^2.*y(1).^2+2.*y(1).*Y+2.*u(3).^2.*x(1).^2+2.*x(1).*X+2.*u(2).^2.*y(1).^2-4.*u(2).^2.*y(1).*Y-u(1).^2.*Y.^2.*u(2).^2+u(2).^2.*u(3).^2.*R.^2-4.*u(1).^2.*x(1).*X+2.*X.*u(1).^4.*x(1)-2.*u(3).^2.*R.^2-u(3).^4.*X.^2-2.*u(3).^2.*X.*u(1).*Y.*u(2)-x(1).^2-y(1).^2+u(3).^4.*R.^2+R.^2+2.*u(1).^2.*x(1).^2-X.^2-Y.^2+4.*x(1).*u(1).*u(2).*y(1)-2.*X.*u(1).^3.*Y.*u(2)+2.*X.*u(1).^3.*u(2).*y(1)+2.*u(1).^3.*x(1).*Y.*u(2)-2.*u(1).^3.*x(1).*u(2).*y(1)+2.*u(1).^2.*Y.*u(2).^2.*y(1)+4.*u(1).*Y.*u(2).*X-4.*u(1).*u(2).*y(1).*X-u(2).^2.*u(1).^2.*x(1).^2+2.*Y.*u(2).^4.*y(1)-2.*u(3).^2.*X.^2.*u(1).^2-2.*u(3).^2.*u(1).^2.*x(1).^2-2.*u(3).^2.*Y.^2.*u(2).^2-2.*u(3).^2.*u(2).^2.*y(1).^2-X.^2.*u(1).^4+2.*X.^2.*u(1).^2-u(1).^4.*x(1).^2-Y.^2.*u(2).^4+2.*Y.^2.*u(2).^2-u(2).^4.*y(1).^2+2.*u(2).^2.*X.*u(1).^2.*x(1)-2.*u(2).^3.*X.*u(1).*Y+2.*u(2).^3.*X.*u(1).*y(1)+2.*u(2).^3.*u(1).*x(1).*Y-2.*u(2).^3.*u(1).*x(1).*y(1)-4.*u(2).*u(1).*x(1).*Y+4.*u(3).^2.*X.*u(1).^2.*x(1)+4.*u(3).^2.*Y.*u(2).^2.*y(1)+2.*u(3).^2.*X.*u(1).*u(2).*y(1)+2.*u(3).^2.*u(1).*x(1).*Y.*u(2)-2.*u(3).^2.*u(1).*x(1).*u(2).*y(1)-u(3).^4.*x(1).^2-u(3).^4.*y(1).^2+2.*u(3).^2.*X.^2+2.*u(3).^2.*Y.^2+u(3).^2.*u(1).^2.*R.^2-u(3).^2.*u(1).^2.*Y.^2-u(2).^2.*u(3).^2.*x(1).^2+2.*u(2).^2.*u(3).^2.*x(1).*X+2.*u(3).^2.*u(1).^2.*y(1).*Y+2.*u(3).^4.*y(1).*Y+2.*u(3).^4.*x(1).*X-u(3).^2.*u(1).^2.*y(1).^2).^(1./2))./(-2.*u(3).^2+u(2).^2.*u(3).^2+1+u(3).^4+u(3).^2.*u(1).^2);

% Filter out imaginary elements (Hemisphere @ P1)
Z1(imag(Z1)~=0) = nan;
% Filter out imaginary elements (Hemisphere @ P2)
Z2(imag(Z2)~=0) = nan;
% Filter out imaginary elements (Cylindar b/w P1,P2)
Zc(imag(Zc)~=0) = nan;

% theta = atan(u(3)/sqrt(u(1)^2 + u(2)^2));
% c = 2*R/tan(theta-pi/2)/2;
% w = u*abs(c)*0;

% % Rotate Before trimming
% theta = atan(u(3)/sqrt(u(1).^2+u(2).^2));
% phi = atan(u(2)/u(1));
% T1 = [cos(phi) -sin(phi) 0;...
%       sin(phi) cos(phi)  0;...
%       0        0         1];
% T2 = [cos(theta)  0 sin(theta);...
%       0           1 0         ;...
%       -sin(theta) 0 cos(theta)];
%   
% T = T1*T2;
% [m,n] = size(X);
% XYZ = [X(:),Y(:),Zc(:)];
% XYZ = (T*XYZ')';
% X = reshape(XYZ(:,1),m,n);
% Y = reshape(XYZ(:,2),m,n);
% Zc = reshape(XYZ(:,3),m,n);

v = [u(2) , -u(1)];
mag = sqrt(v(1)^2 + v(2)^2);
mag = max(eps,mag);
v = R*v/mag;
boxX = [-v(1) + x(1), v(1) + x(1), v(1) + x(2), -v(1) + x(2)]; 
boxY = [-v(2) + y(1), v(2) + y(1), v(2) + y(2), -v(2) + y(2)]; 
% line(boxX,boxY,'color','r')
in = inpolygon(X,Y,boxX,boxY);


% xy = sqrt(x.^2 + y.^2);
% w = [u(3) -sqrt(u(1)^2 + u(2)^2)];
% mag = sqrt(w(1)^2 + w(2)^2);
% mag = max(eps,mag);
% w = R*w/mag;
% 
% boxXY = [xy(1) + w(1) xy(1) xy(2) xy(2) + w(1)];
% boxZ  = [ z(1) + w(2)  z(1)  z(2)  z(2) + w(2)];
% in2 = inpolygon(sqrt(X.^2 +Y.^2),Zc,boxXY,boxZ);
% figure(2345)
% line(boxXY,zeros(size(boxZ))-1,boxZ,'color','r')

Zc(~in) = nan;

% % Rotate Back
% XYZ = [X(:),Y(:),Zc(:)];
% XYZ = (inv(T)*XYZ')';
% X = reshape(XYZ(:,1),m,n);
% Y = reshape(XYZ(:,2),m,n);
% Zc = reshape(XYZ(:,3),m,n);

% Some boolean operations
Z = min(Z1,Z2);
Z = min(Z,Zc);




function createBottom(minX,maxX,minY,maxY,minZ)

X = zeros(4);
Y = X;
Z = X;

X(:,1:2) = minX;
X(:,3:4) = maxX;
Y(1:2,:) = minY;
Y(3:4,:) = maxY;
Z(2:3,2:3) = minZ;

hold on
h = surf(X,Y,Z,...
    'edgecolor','none',...
    'facecolor',[.2 .8 1],...
    'facelight','gouraud',...
    'specularstrength',0,...
    'ambientstrength',.3);
