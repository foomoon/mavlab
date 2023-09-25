function [x,y,z,ind] = surf2zigzag(X,Y,Z,step,res)
% ZIGZAG - Creates a zigzag line on a surface
%   [x,y,z] = zigzag(X,Y,Z,step,res) inputs surface data (X,Y,Z), stepover
%   distance between parallel cuts (step), and resolution (res). Resolution
%   is the inverse of the number of data points per unit of measurement.
%

global dcDebug

PtsPerInch = 1/res;

% Get Leading/Trailing Edge from Surface Data
[e1,e2] = surf2edge(X,Y,Z);
x1 = e1(:,1); y1 = e1(:,2); z1 = e1(:,3);
x2 = e2(:,1); y2 = e2(:,2); z2 = e2(:,3);

% Grid from min(x) to max(x)
xext  = max(x2)-min(x2);
% xnpts = round(xext/(step) + 1);
xnpts0 = length(x2);
xii = linspace(min(x2),max(x2),xnpts0);

% Interpolate leading/trailing edges onto grid
method = 'linear';
yl = interp1(x1,y1,xii,method); % Leading  Edge
yt = interp1(x2,y2,xii,method); % Trailing Edge
zl = interp1(x1,z1,xii,method); % Leading  Edge
zt = interp1(x2,z2,xii,method); % Trailing Edge
% figure; line(xi,yl,zl,'color','g','marker','.');line(xi,yt,zt,'color','r','marker','.'); axis equal
isZeroL = (yl-yt == 0);
xii(isZeroL) = [];
yl(isZeroL) = [];
zl(isZeroL) = [];
yt(isZeroL) = [];
zt(isZeroL) = [];
xnpts = round(xext/(step) + 1);
xi = linspace(min(x2),max(x2),xnpts);
yl = csapi(xii,yl,xi); % Leading  Edge
yt = csapi(xii,yt,xi); % Trailing Edge
zl = csapi(xii,zl,xi); % Leading  Edge
zt = csapi(xii,zt,xi); % Trailing Edge
% figure; line(xi2,yl,zl,'color','g','marker','.');line(xi2,yt,zt,'color','r','marker','.'); axis equal
% Remove any points where the distance between leading and trailing edge is
% zero

% isZeroL = (yl-yt == 0);
% xi(isZeroL) = [];
% yl(isZeroL) = [];
% zl(isZeroL) = [];
% yt(isZeroL) = [];
% zt(isZeroL) = [];

maxiter = length(xi);
x = []; y = []; zfix = zeros(maxiter,2); ind = zfix;
for i=1:maxiter
    yext  = yl(i)-yt(i);
    ynpts = yext*PtsPerInch;
    if mod(i,2)
        ytemp = linspace(yl(i),yt(i),ynpts)';
        zfix(i,:) = [zl(i) zt(i)];
    else
        ytemp = linspace(yt(i),yl(i),ynpts)';
        zfix(i,:) = [zt(i) zl(i)];
    end    
    xtemp = xi(i)*ones(size(ytemp));
    ind(i,:) = [length(x)+1 length(x)+length(xtemp)];
    x = [x; xtemp];
    y = [y; ytemp];    
    workbar(i/maxiter,'Computing ZigZag Path');
end

% interpolate z data
z = griddata(X,Y,Z,x,y,'cubic');
% z = interp2(X,Y,Z,x,y,'cubic');
z(ind(:)) = zfix(:);
ind = ind(:);

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  Notice: %0.0f fixed toolpath point(s)\n',length(ind))
end


