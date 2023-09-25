function nrbOffset(srf,res,r)

if nargin < 3,
    r = 1;
end

% Parametric coordinates s,t
s = linspace(0,1,res(1));
t = linspace(0,1,res(2));


% Evaluate position of actual surface at s,t
p = nrbeval(srf,{s t});
x = squeeze(p(1,:,:));
y = squeeze(p(2,:,:));
z = squeeze(p(3,:,:));

% NURBS Derivative
dsrf = nrbderiv(srf);

% Evaluate derivative at s,t
[p1, dp] = nrbdeval(srf, dsrf, {s, t});

% Normalize Tangent/Binormal
T = vecnorm(dp{1});
B = vecnorm(dp{2});
% Compute Normal via cross product and normalize it
N = veccross(T,B);
N = vecnorm(N);
% Scale Normal to r
N = N*r;

% figure; 
if any(size(x))==1,
    line(x,y,z)
else
    surf(x,y,z)
end
hold on
plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
quiver3(p1(1,:),p1(2,:),p1(3,:),N(1,:),N(2,:),N(3,:),0);
axis equal
hold off

