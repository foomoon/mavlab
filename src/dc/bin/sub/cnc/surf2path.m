function [x,y,z] = surf2path(X,Y,Z,D,step,tol,type,flip,dir),
% SURF2PATH Wrapper for surf2spiral and surf2zigzag
%
%   [x,y,z] = surf2path(X,Y,Z,D,step,tol,type)
%
%   See also surf2spiral surf2zigzag

global dcDebug

warning off

% Calculate resolution of point spacing
Rmin = minCurvature(X,Y,Z);
res = tol; %tol2res(tol,Rmin);

% Create offset surface
if flip
    [Xo,Yo,Zo,v] = offsetSurface(X',Y',Z',D/2);
    Xo = Xo'; Yo = Yo'; Zo = Zo'; v = v';
%     disp('Surf2path: Flipping')
else
    [Xo,Yo,Zo,v] = offsetSurface(X,Y,Z,D/2);
end

% Toolpath selection
switch type,
    case 'spiral',
        [x,y,z] = surf2spiral(Xo,Yo,Zo,step,res);
    case 'zigzag',
        switch dir(1)
            case 'h'
                [y,x,z,i] = surf2zigzag(Yo,Xo,Zo,step,res);
            case 'v'
                [x,y,z,i] = surf2zigzag(Xo,Yo,Zo,step,res);
            otherwise
                [x,y,z,i] = surf2zigzag(Xo,Yo,Zo,step,res);
        end
%         Zo(v) = nan;
%         z2 = griddata(Xo,Yo,Zo,x,y,'cubic');
%         z2(i) = 0;
%         z2 = z2-z2;
%         z = z + z2;
    otherwise
end

warning on

[x,y,z] = prepToolpath(x,y,z,step,res);

if dcDebug
    L = sum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));
    n = length(x);
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  %s toolpath %0.0f point(s)\n',type,n)
    fprintf(1,'  %s toolpath %0.1f long\n',type,L)
end
