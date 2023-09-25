function [X,Y,Z] = nrbevalany(nrb,npts)

% [m,n] = size(x);

% convert NURBS to B-spline
nrb.form = 'B-';
nrb.coefs(4,:,:) = [];
nrb.dim = 3;


% convert form to piece-wise polynomial
pp = fn2fm(nrb,'pp');


% Evaluate on interval
% xyz = ppual(pp,{x(:) y(:)});
% z=ppual(pp,[x(:) y(:)]');
% z = xyz(3,:);
% z = reshape(z,m,n);

b = pp.breaks;

if nargin < 2
    npts = [51 51];
else
    npts = npts + 1;
end

ext = .05;
start = 0 - ext;
finish = 1 + ext;
xx = {linspace(start,finish,npts(1)); linspace(start,finish,npts(2))};
% xx = {linspace(b{1}(1),b{1}(end),npts(1)) linspace(b{2}(1),b{2}(end),npts(2))};
% xx = [linspace(b{1}(1),b{1}(end),npts); linspace(b{2}(1),b{2}(end),npts)];
v = ppual(pp,xx);
v = permute(reshape(v,[3,npts(1),npts(2)]),[2 3 1]);
% v = permute(reshape(v,[3,length(xx(1,:)),length(xx(2,:))]),[2 3 1]);
X = v(:,:,1);
Y = v(:,:,2);
Z = v(:,:,3);

% Debug
surf(X,Y,Z)

