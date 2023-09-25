function out = gnurbs(mode,res)
% GNURBS  Interactive manipulation of NURBS surface/curve
%
%    For use with the NURBS toolbox, GNURBS allows the user to have an
%    intuitive manipulation of a NURBS surface/curve via control points.  
%
%    GNURBS with no inputs will generate a test surface for manipulation.
%    GNURBS(NRB), where NRB is a structure defining a NURBS surface/curve
%    (as defined in the NURBS toolbox).
%    GNURBS(NRB,RES), where RES is a 2 element vector defines the
%    resolution of surface/curve rendering in the u and v directions, 
%    ie RES = [URES,VRES].  The default is [20,20] if no argument is
%    specified.
%    
%    <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=312&objectType=file">NURBS toolbox</a> can be downloaded for free at
%    http://www.mathworks.com/matlabcentral/fileexchange/
%
%    Example usage:
%    % Generate a NURBS surface
%    gnurbs



if ~nargin
    mode = 'start';
    nrb = nrbtestsrf; % Default Surface
    res = [20 20]; % Default Resolution
elseif isstruct(mode)
    nrb=mode;
    mode='start';
    if nargin<2
        res = [20 20]; % Default Resolution
    end
elseif isempty(mode)
    nrb=nrbtestsrf;
    mode='start';
    if nargin<2
        res = [20 20]; % Default Resolution
    end
elseif ischar(mode)
    fig = gcf;
    ud = get(gcbo,'userdata');
%     disp(mode)
else
    error('GNURBS: Wrong type of input')
end


switch mode
    case 'start'
        % Build surface from NURBS structure
        [x,y,z] = nrb2xyz(nrb,res);
        % Get control net
        [xn,yn,zn] = nrb2net(nrb);
        % Plot
        fig = gcf;
        ax = gca;
        axis equal
        
        n = length(nrb.order);
        
        if n > 1
            hsurf = surface(x,y,z,'parent',ax,'facecolor','interp',...
                        'facealpha',.5,'edgealpha',.2);
            hnet  = surface(xn,yn,zn,'parent',ax,'facecolor','none',...
                        'edgecolor','k','marker','o','linestyle',':',...
                        'markeredgecolor','k','markerfacecolor','r',...
                        'markersize',5);
        else            
            hsurf = line(x,y,z,'parent',ax,'linestyle','-','color','k');
            hnet  = line(xn,yn,zn,'parent',ax,...
                        'color','k','marker','o','linestyle',':',...
                        'markeredgecolor','k','markerfacecolor','r',...
                        'markersize',5);                       
        end            
        hcur  = line(0,0,0,'marker','o','markerfacecolor','g',...
                     'markeredgecolor','k','hittest','off','visible','off',...
                     'markersize',5);
                    

      
        % Store handles and data in figure userdata
        ud.ax    = ax;
        ud.hsurf = hsurf;
        ud.hnet  = hnet;
        ud.hcur  = hcur;
        ud.nrb   = nrb;
        ud.res   = res;
        set(hnet,'UserData',ud);
        
        % Set Callbacks
        set(hsurf,'ButtonDownFcn','','hittest','off');
        set(hnet, 'ButtonDownFcn','gnurbs(''click'')');
        
        % Outputs
        if nargout
            out = hsurf;
        end
        
    case 'click'        
        ud.ind = getClosestPt(ud.hnet);
        [xn,yn,zn] = nrb2net(ud.nrb);
        setxyz(ud.hcur,xn(ud.ind),yn(ud.ind),zn(ud.ind));
        set(ud.hcur,'visible','on');  
        storeparams(fig,ud);
        set(fig,'WindowButtonMotionFcn','gnurbs(''move'')');
        set(fig,'WindowButtonUpFcn','gnurbs(''up'')');
%         set(fig,'UserData',ud);
        
    case 'move'
        currpt = get(ud.ax,'CurrentPoint');
        [xn,yn,zn] = nrb2net(ud.nrb);
        [x,y,z] = viewptx(currpt,[xn(ud.ind),yn(ud.ind),zn(ud.ind)]);
        xn(ud.ind) = x;
        yn(ud.ind) = y;
        zn(ud.ind) = z;
        coefs = net2nrb(xn,yn,zn);
        ud.nrb.coefs = coefs;
        [X,Y,Z] = nrb2xyz(ud.nrb,ud.res);
        setxyz(ud.hcur,x,y,z);     % Update Hightlight point
        setxyz(ud.hsurf,X,Y,Z);    % Update Surface
        setxyz(ud.hnet,xn,yn,zn);  % Update Control Net
        setxyz(ud.hcur,x,y,z);     % Update Hightlight point
        set(fig,'UserData',ud);
        set(ud.hnet,'userdata',ud);
        set(ud.hsurf,'userdata',ud.nrb);
    case 'up'
        set(ud.hcur,'visible','off');
        resetparams(fig);
    otherwise
        error('GNURBS: Wrong input string')
end






function [X,Y,Z] = viewptx(vdir,pt)
% VIEWPTX  Finds the intersection point of mouse click and viewing plane
%   The viewing plane is defined as the plane normal to the current viewing
%   angle in the current figure.  When the current view changes using
%   Rotate3d or View(az,el), Viewptx adjusts the viewing plane accordingly.
%   Since the position of the viewing plane is ambiguous, it is necessary
%   to define a second argument that specifies a point the plane passes 
%   through in euclidean space.
%
%   [X,Y,Z] = VIEWPTX(VDIR,PT) where VDIR is a 3x2 array specifying the
%   endpoints of a line segment, or "stabbing vector".  This value can be
%   attained by VDIR = GET(GCA,'CURRENTPOINT').  The second arguement, PT,
%   is a point the plane passes through.  
%
%   VIEWPTX is a helper function, not intended to be used on its own.
%   However, it could be adapted for other uses.  In it's current state, it
%   is designed to allow a user to graphically manipulate 3D points in a
%   useful and intuitive way.  The end behavior of this function is to
%   translate data in the plane normal to the viewing angle and passing
%   through a specified point.
%
%   See also: VIEW

% Daniel Claxton
% 02-Mar-2007
% dclaxton@ufl.edu


M = view;                          % Current View matrix
R = inv(M(1:3,1:3));               % Create Rotation matrix from current view
up = [[1 0 0]; [0 1 0]; [0 0 0]];  % Arbitrary Data points on plane with normal [0 0 1]
T = R*up;                          % Perform Rotation of "up" data points
x = T(1,:) + pt(1);                % Translate X
y = T(2,:) + pt(2);                % Translate Y
z = T(3,:) + pt(3);                % Translate Z

% Concatenate vdir onto x,y,z to make expressions easier to read
x = [x vdir(1,1) vdir(2,1)];       
y = [y vdir(1,2) vdir(2,2)];
z = [z vdir(1,3) vdir(2,3)];

% Fun Numerical Stuff (Solve System)
num = det([ones(1,4); x(1:4); y(1:4); z(1:4)]);
den = det([[ones(1,3); x(1:3); y(1:3); z(1:3)] [0; (x(5)-x(4)); (y(5)-y(4)); (z(5)-z(4))]]);

t = -num/den;

X = x(4) + (x(5)-x(4))*t ;
Y = y(4) + (y(5)-y(4))*t ;
Z = z(4) + (z(5)-z(4))*t ;




function ind = getClosestPt(h)
% GETCLOSESTPT  Finds the closest point in a data set to a line
%    In this case, the data set is defined by a surface and the line is
%    simply the result of get(gca,'CurrentPoint').
%    I = GETCLOSESTPT(H) returns the index of the point in H closest to a 
%    mouse click in the figure.  H is a handle to a surface object.

% Daniel Claxton
% 02-Mar-2007
% dclaxton@ufl.edu

x = get(h,'xdata');
y = get(h,'ydata');
z = get(h,'zdata');
if isempty(z), z = x*0; end
n = numel(x);
ax = get(h,'parent');
cpt = get(ax,'CurrentPoint');
% Find vector in direction of line
u = diff(cpt);
% Normalize line vector
mag = sqrt(sum(u.^2));
mag = max(eps,mag);
u = u/mag;
% First Point on line (P1)
v = cpt(1,:);
d = zeros(n,1);
for i=1:n
    pt = [x(i) y(i) z(i)];
    % Make vector between P1 and point
    w = pt - v;
    % Dot product gives us the magnitude of projection of w onto line
    f = u*w';
    % Multiply magnitude of projection times unit vector u (direction of line)
    % Shift from orgin to P1
    p = f*u + v;
    % Find shortest distance from point to point on line
    d(i) = sqrt(sum((pt-p).^2,2));
end
[dist,ind] = min(d);


function [x,y,z] = nrb2net(srf)
% NRB2NET  Get control net from NURBS structure
x=squeeze(srf.coefs(1,:,:));
y=squeeze(srf.coefs(2,:,:));
z=squeeze(srf.coefs(3,:,:));


function nrb = net2nrb(x,y,z)
% NET2NRB Get NURBS coefficients from control net
[m,n]=size(x);
x = reshape(x,1,m,n);
y = reshape(y,1,m,n);
z = reshape(z,1,m,n);
one = ones(1,m,n);

nrb = cat(1,x,y);
nrb = cat(1,nrb,z);
nrb = cat(1,nrb,one);


function [x,y,z] = nrb2xyz(nrb,res)
% NRB2XYZ  Evaluate NURBS with a given resolution

n = length(nrb.order);
if n > 1
    nu = res(1);
    nv = res(2);
    p = nrbeval(nrb,{linspace(0,1,nu),linspace(0,1,nv)});
else
    p = nrbeval(nrb,linspace(0,1,res(1)));
end
x = squeeze(p(1,:,:));
y = squeeze(p(2,:,:));
z = squeeze(p(3,:,:));


function setxyz(h,x,y,z)
% SETXYZ  Set X,Y,Z Data of graphics handle object
set(h,'Xdata',x);
set(h,'Ydata',y);
set(h,'Zdata',z);

if strcmp(get(h,'type'),'surface')
    set(h,'cdata',z/max(z(:)));
end
                  

function ud = storeparams(fig,ud)
ud.udata = get(fig,'userdata');
ud.down  = get(fig,'windowbuttondownfcn');
ud.move  = get(fig,'windowbuttonmotionfcn');
ud.up    = get(fig,'windowbuttonupfcn');
set(fig,'userdata',ud);

function resetparams(fig)
ud = get(fig,'userdata');
set(fig,'userdata',ud.udata','windowbuttondownfcn',ud.down,...
    'windowbuttonmotionfcn',ud.move,'windowbuttonupfcn',ud.up);