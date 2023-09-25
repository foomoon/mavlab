function [x,y,z] = surf2spiral(X,Y,Z,step,res)

% Limit res
res0 = res;
res = min(res,.5*step);

% Get Border Data 
x1 = fliplr(X(1,:)); x2 =  X(:,1)'; x3 = X(end,:); x4 = flipud(X(:,end))';
y1 = fliplr(Y(1,:)); y2 =  Y(:,1)'; y3 = Y(end,:); y4 = flipud(Y(:,end))';
z1 = fliplr(Z(1,:)); z2 =  Z(:,1)'; z3 = Z(end,:); z4 = flipud(Z(:,end))';

% Make sure normals always point "in"
i=round(length(x2)/2); 
dir = (x2(i) - x2(i-1)) * (y2(i+1) - y2(i)) - (y2(i) - y2(i-1)) * (x2(i+1) - x2(i));
dir = -dir/abs(dir);
step = step*dir;

sp1 = splineval(x1,y1,z1,res); % Left   Edge (D)
sp2 = splineval(x2,y2,z2,res); % Top    Edge (A)
sp3 = splineval(x3,y3,z3,res); % Right  Edge (B)
sp4 = splineval(x4,y4,z4,res); % Bottom Edge (C)

% xext = max(X(:))-min(X(:));
% yext = max(Y(:))-min(Y(:));
xext = max(X(:))-min(X(:));
xi = linspace(min(X(:)),max(X(:)),round(1/res));
yl = interp1(sp2(:,1),sp2(:,2),xi,'cubic');
yt = interp1(sp4(:,1),sp4(:,2),xi,'cubic');
yext = abs(max(yl-yt));

xscale = 2*step/xext;
yscale = 2*step/yext;

maxXiter = 1/(xscale);
maxYiter = 1/(yscale);
maxiter = floor(min(maxXiter,maxYiter));

A = [];
B = [];
D = [];

sp1(:,3) = []; sp2(:,3) = []; sp3(:,3) = []; sp4(:,3) = [];

warning off
for i=1:maxiter
    % Find normals to previous curve;
    sp1 = step*cnorm(sp1) + sp1;
    sp2 = step*cnorm(sp2) + sp2;
    sp3 = step*cnorm(sp3) + sp3;
    sp4 = step*cnorm(sp4) + sp4;
    
    [pt,ind] = xsect(sp2,sp4); % A n C
    if ~isempty(ind)
        % Won't Occur Often... but will likely occur in last few iterations
        sp2 = trimline(sp2,pt,ind(:,1));
        sp4 = trimline(sp4,pt,ind(:,2));
        sp1 = []; sp3 = [];
    else
        [pt,ind] = xsect(sp1,sp3); % B n D
        if ~isempty(ind)
            % Won't occur often...
            sp1 = trimline(sp1,pt,ind(:,1));
            sp3 = trimline(sp3,pt,ind(:,2));
            sp2 = []; sp4 = [];
        else
            % If no B n D
            % If no A n C
            % Will Likely occur most of the time
            [sp1,sp2] = trimlines(sp1,sp2); % D n A
            [sp2,sp3] = trimlines(sp2,sp3); % A n B
            [sp3,sp4] = trimlines(sp3,sp4); % B n C
            [sp1,sp4] = trimlines(sp1,sp4); % C n D
        end

    end
    
        % Re-Spline
    sp11 = splineval(sp1,res0); % Left   Edge (D)
    sp22 = splineval(sp2,res0); % Top    Edge (A)
    sp33 = splineval(sp3,res0); % Right  Edge (B)
    sp44 = splineval(sp4,res0); % Bottom Edge (C)
    
    % Re-Spline
    sp1 = splineval(sp1,res); % Left   Edge (D)
    sp2 = splineval(sp2,res); % Top    Edge (A)
    sp3 = splineval(sp3,res); % Right  Edge (B)
    sp4 = splineval(sp4,res); % Bottom Edge (C)

    A = [sp1; sp2; sp3; sp4];    
    
    D = [sp11; sp22; sp33; sp44];
    B = [B; D];


   workbar(i/maxiter,'Creating Spiral Path'); 
end
warning on


sp1 = splineval(x1,y1,z1,res0); % Left   Edge (D)
sp2 = splineval(x2,y2,z2,res0); % Top    Edge (A)
sp3 = splineval(x3,y3,z3,res0); % Right  Edge (B)
sp4 = splineval(x4,y4,z4,res0); % Bottom Edge (C)
C = [sp1; sp2; sp3; sp4];

B(:,3) = griddata(X,Y,Z,B(:,1),B(:,2),'cubic');
C = [C; B];
x = C(:,1); y = C(:,2); z = C(:,3);




function [sp1,sp2] = trimlines(sp1,sp2)
[pt,ind] = xsect(sp1,sp2); % D n A
if isempty(ind), return, end
sp1 = trimline(sp1,pt,ind(:,1));
sp2 = trimline(sp2,pt,ind(:,2));

function N = cnorm(s)

if isempty(s)
    N = [];
    return
end
x = s(:,1);
y = s(:,2);

dx = gradient(x);
dy = gradient(y);

N = [dy -dx];
Nmag = sqrt(sum(N.^2,2));
N = N./Nmag(:,[1 1]);





function out = splineval(x,y,z,res)
if isempty(x) | isnan(sum(x(:)))
    out = [];
    return
end
if nargin == 2,
    res = y;
    y = x(:,2);
    x = x(:,1);
    if length(x) >= 4
        [x,y]=rmselfx(x,y);
    end
    z = zeros(size(x));
end

rm = find(isnan(x(:)+y(:)+z(:)));
x(rm) = [];
y(rm) = [];
z(rm) = [];

% % remove repeated points
% rmrep = find((diff(sqrt(x.^2 + y.^2 + z.^2)))==0);
% x(rmrep) = [];
% y(rmrep) = [];
% z(rmrep) = [];
% if length(x) < 2,
%     out = [x(:) y(:) z(:)];
% else
ctrlpts = [x(:)'; y(:)'; z(:)'];
L = sum(sqrt(sum(diff(ctrlpts,1,2).^2)));
npts = round(L/res);
sp = cscvn(ctrlpts);
last = sp.breaks(end);
val = linspace(0,last,npts);
out = fnval(sp,val)';
% end
if nargin == 2 & ~isempty(out)
    out = out(:,1:2);
end

function [pt,ind] = xsect(p1,p2)
if size(p1,1) < 1 | size(p2,1) < 1
    pt = [];
    ind = [];
    return
end

[x0,y0] = curveIntersect(p1(:,1),p1(:,2),p2(:,1),p2(:,2));
ind = [];
for i=1:length(x0)
    d1 = sqrt( (p1(:,1)-x0(i)).^2 + (p1(:,2)-y0(i)).^2 );
    d2 = sqrt( (p2(:,1)-x0(i)).^2 + (p2(:,2)-y0(i)).^2 );
    [m,ind(i,1)] = min(d1);
    [m,ind(i,2)] = min(d2);
end
pt = [x0(:) y0(:)];


function p = trimline(p,pt,ind)

p = p(:,1:2);
L = lineLength(p);
for i=1:length(ind)
    pA = [pt(i,:); p(ind(i)+1:end,:)];
    pB = [p(1:ind(i)-1,:); pt(i,:)];
    L1 = lineLength(pA)/L;
    L2 = (L-lineLength(pB))/L;
    L1 = max(L1,L2);
    if L1 > .5
        p(1:ind(i),:) = nan;
        p(ind(i),:) = pt(i,:);
    else
        p(ind(i):end,:) = nan;
        p(ind(i),:) = pt(i,:);
    end
end
p(isnan(p(:,1)),:) = [];


function L = lineLength(p)
L = sum(sqrt(sum(diff(p).^2,2)));