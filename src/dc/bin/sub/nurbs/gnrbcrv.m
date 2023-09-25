function h = gnrbcrv(fig)
% gnrbcrv - Interactively draw interpolating spline similar to gline
%   gnrbcrv(FIG) draws a spline by left - clicking the mouse at the control
%   points of the spline in the figure FIG. Right clicking will end the
%   spline.
% 
%   H = gnrbcrv(FIG) Returns the handle to the line.
% 
%   gnrbcrv with no input arguments draws in the current figure.

% erasemode = 'xor';
erasemode = 'normal';

if nargin<1, 
  draw_fig = gcf;
  fig = draw_fig;
  s = 'start'; 
  gca;
  ax = get(draw_fig,'CurrentAxes');
%   set(draw_fig,'color','k');
%   set(ax,'color','k');
%   op = get(ax,'outerposition');
%   set(ax,'position',op);
%   set(draw_fig,'toolbar','none','menu','none');
end


if isstr(fig), 
   s = fig;
   draw_fig = gcbf;
   ud = get(draw_fig,'UserData');
else
   isline = strcmp('line',get(fig,'type'));
   s = 'start';
   if isline       
       draw_fig = gcf;
   else
       draw_fig = fig;
   end
end

ax = get(draw_fig,'CurrentAxes');
if isempty(ax)
   ax = axes('Parent',draw_fig);
end

gcacolor = get(ax,'Color');

% Check to see if we're clicking to exit
mouse = lower(get(draw_fig,'selectiontype'));
if strcmp(mouse(1:3),'alt') & strcmp(s,'down')
    s = 'last';
end
% Normal:    Click left mouse button.
% Extend:    Shift - click left mouse button or click both left and right 
%            mouse buttons,  or click middle mouse button. 
% Alternate: Control - click left mouse button or click right mouse button.
% Open:      Double-click any mouse button.


switch s
   case 'start'
   oldtag = get(draw_fig,'Tag');
   figure(draw_fig);
   if any(get(ax,'view')~=[0 90]), 
     set(draw_fig, 'Tag', oldtag);
     error('stats:gnrbcrv:NotTwoDimensional','gnrbcrv works only for 2-D plots.');
   end
   
   xlimits = get(ax,'Xlim');
   x = (xlimits + flipud(xlimits))./2;
   ylimits = get(ax,'Ylim');
   y = (ylimits + flipud(ylimits))./2;
   
   if ~isline
       % Initialize spline
       hline = line(x,y,'Parent',ax,'Visible','off','eraseMode',erasemode,...
           'color',1-get(ax,'color'),'clipping','off','hittest','off'); 
       npts = 1;
   else
       hline = fig;
       for i=1:length(hline)
           x = get(hline(i),'xdata');
           y = get(hline(i),'ydata');
           z = get(hline(i),'zdata');
           if isempty(z), z = x*0; end
           [xspline,yspline,zspline] = splineInterp(x,y,z);
           set(hline(i),'userdata',[x(:) y(:)],'xdata',xspline(:),'ydata',yspline(:),'zdata',zspline(:));
           set(hline(i),'buttondownfcn','gnrbcrv(''edit'')')
       end
       return
   end
   
   % Change pointer to cross hair
   set(draw_fig,'Pointer','cross');
   
   % Intialize Current location text box
   Pt = get(ax,'CurrentPoint');   
   Pt = Pt(1,1:2);      
   htext = text(0,0,' ');
   pos = get(htext,'position');   
   pos(1:2) = Pt + 0.03*[abs(diff(xlimits)) -abs(diff(ylimits))];
   curloctext = sprintf(' X: %1.4f\n Y: %1.4f ',Pt);   
   set(htext,'string',curloctext,'position',pos)
   set(htext,'fontsize',8,'edgecolor',[.5 .5 .5],'Backgroundcolor',[.98,.98,.98])

   % Save current window functions and data
   bdown = get(draw_fig,'WindowButtonDownFcn');
   bup = get(draw_fig,'WindowButtonUpFcn');
   bmotion = get(draw_fig,'WindowButtonMotionFcn');
   oldud = get(draw_fig,'UserData');
   
   % Create new window functions
   set(draw_fig,'WindowButtonUpFcn','gnrbcrv(''first'')')
   set(draw_fig,'WindowButtonMotionFcn','gnrbcrv(''motion'')')
   set(draw_fig,'WindowButtonDownFcn','')
   
   set(draw_fig,'doublebuffer','on')
 
   % Save drawing data as in 'UserData'
   ud.hline = hline;
   ud.htext = htext;
   ud.pts = [];
   ud.buttonfcn = {bdown; bup; bmotion};
   ud.oldud = oldud;
   ud.oldtag = oldtag;
   ud.npts = npts;
   ud.xlimits = xlimits;
   ud.ylimits = ylimits;
   set(draw_fig,'UserData',ud);
   
   if nargout == 1
      h = hline;
   end

case 'motion'
   
   Pt2 = getCurPt(ax);  
   pos = get(ud.htext,'position');
   
   pos(1:2) = Pt2 + [0.03 0.07].*[abs(diff(ud.xlimits)) -abs(diff(ud.ylimits))];
   curloctext = sprintf(' X: %1.4f\n Y: %1.4f ',Pt2);
   set(ud.htext,'string',curloctext,'position',pos)
   
  
   if isempty(ud.pts);
      return;
   end

   [xspline,yspline,zspline] = splineInterp(ud.pts(:,1),ud.pts(:,2));   
   set(ud.hline,'Xdata',xspline,'Ydata',yspline,'Zdata',zspline, ...
        'linestyle','-', 'linewidth', 1, 'Color',1-gcacolor,'Visible','on'); 
   ud.pts(ud.npts+1,:) = Pt2;
   set(draw_fig,'UserData',ud);

case 'first'   
   Pt1 = get(ax,'CurrentPoint'); 
   ud.pts = [Pt1(1,1:2); Pt1(1,1:2)];

   set(draw_fig,'WindowButtonDownFcn','gnrbcrv(''down'')','UserData',ud)
   
case 'down'         
   Pt1 = get(ax,'CurrentPoint'); 
   ud.pts = [ud.pts; Pt1(1,1:2)]; 
   ud.npts = ud.npts + 1;   

   set(draw_fig,'WindowButtonUpFcn','');
   set(draw_fig,'WindowButtonDownFcn','gnrbcrv(''down'')','UserData',ud)
      
case 'last'
   delete(ud.htext);
   bfcns = ud.buttonfcn;
   set(draw_fig,'windowbuttondownfcn',bfcns{1},'windowbuttonupfcn',bfcns{2}, ...
         'windowbuttonmotionfcn',bfcns{3},'Pointer','arrow', ...
		 'Tag',ud.oldtag,'UserData',ud.oldud)
   set(ud.hline,'UserData',ud.pts,'ButtonDownFcn','gnrbcrv(''edit'')','hittest','on')

case 'edit'
   
   % Get anchor points
   hline = gcbo;   
   pts = get(hline,'userdata');
   hmark = line(pts(:,1),pts(:,2));
   
%    set(hline,'hittest','off')
   set(hmark,'linestyle',':','marker','s','markersize',4,'clipping','off')
   set(hmark,'markeredgecolor','b','markerfacecolor','k','handlevisibility','off')
   
   xlimits = get(ax,'Xlim');
   ylimits = get(ax,'Ylim');
  
    % Save current window functions and data
   bdown = get(draw_fig,'WindowButtonDownFcn');
   bup = get(draw_fig,'WindowButtonUpFcn');
   bmotion = get(draw_fig,'WindowButtonMotionFcn');
   oldud = get(draw_fig,'UserData');
   pointer = get(draw_fig,'Pointer');
   units = get(draw_fig,'Units');
   
   ud.hline = hline;
   ud.hmark = hmark;
   ud.oldud = oldud;
   ud.bdown = bdown;
   ud.bup = bup;
   ud.bmotion = bmotion;
   ud.pointer = pointer;
   ud.units = units;
   ud.xlimits = xlimits;
   ud.ylimits = ylimits;
   set(draw_fig,'userdata',ud)
   
   % Set defaults
   set(draw_fig,'WindowButtonDownFcn','')
   set(draw_fig,'WindowButtonUpFcn','')
   set(draw_fig,'WindowButtonMotionFcn','')
   set(draw_fig,'pointer','hand')
   set(draw_fig,'units','normalized')
   
   set(ud.hmark,'ButtonDownFcn','gnrbcrv(''edit2'')')
   set(ud.hline,'ButtonDownFcn','gnrbcrv(''editinsert'')')
   
case 'editinsert'
    pts = get(ud.hline,'userdata');
    [m,n] = size(pts);    if m<n, pts = transpose(pts); end
    pt = get(ax,'CurrentPoint');
    pt = pt(1,1:2);
    pts = insertPoint(pts,pt);
    set(ud.hmark,'xdata',pts(:,1),'ydata',pts(:,2))
    [xspline,yspline,zspline] = splineInterp(pts(:,1),pts(:,2));
    set(ud.hline,'userdata',pts,'xdata',xspline,'ydata',yspline,'zdata',zspline)
   

   
case 'edit2'
   % Get control point closest to mouse click
%    Pt = get(ax,'CurrentPoint'); 
%    Ptx = Pt(1,1);
%    Pty = Pt(1,2);
%    x = get(ud.hmark,'xdata');
%    y = get(ud.hmark,'ydata');
%    d = sqrt( (Ptx-x).^2 + (Pty-y).^2 );
%    [m,i] = min(d);
   ud.ind = getClosestPt(ud.hmark);
   if strcmp(mouse(1:3),'alt'),
%       disp('deleting point')
      x(i)=[];
      y(i)=[];
      [xspline,yspline,zspline] = splineInterp(x,y);
      set(ud.hmark,'xdata',x,'ydata',y);
      set(ud.hline,'xdata',xspline,'ydata',yspline,'zdata',zspline,'userdata',[x(:) y(:)])
      return
   end
   set(draw_fig,'windowbuttonmotionfcn','gnrbcrv(''editmove'')',...
       'windowbuttonupfcn','gnrbcrv(''editmoveend'')','userdata',ud)%
%    disp('Edit point')
   
case 'editmove'
    pts = get(ax,'currentpoint');
    x = get(ud.hmark,'xdata'); 
    y = get(ud.hmark,'ydata'); 
    z = get(ud.hmark,'zdata');
    
    xx = x(ud.ind); yy = y(ud.ind);
    if isempty(z), z=x*0; end; zz = z(ud.ind);
    [X,Y,Z] = plnlinex(pts,xx,yy,zz);
%    Pt = getCurPt(ax);

   x(ud.ind) = X;
   y(ud.ind) = Y;    
   z(ud.ind) = Z;
   set(ud.hmark,'xdata',x,'ydata',y,'zdata',z);
   [xspline,yspline,zspline] = splineInterp(x,y,z);
   set(ud.hline,'xdata',xspline,'ydata',yspline,'zdata',zspline,'userdata',[x(:) y(:) z(:)])
   
case 'editmoveend'    
    set(draw_fig,'windowbuttonmotionfcn','')
    set(ax,'buttondownfcn','gnrbcrv(''editmoveend_final'')')
case 'editmoveend_final'
    % Return figure to initial settings
    set(ax,'buttondownfcn','')
    set(draw_fig,'windowbuttonmotionfcn',ud.bmotion,'windowbuttonupfcn',ud.bup,...
        'windowbuttondownfcn',ud.bdown,'userdata',ud.oldud,'pointer',ud.pointer,...
        'units',ud.units) 
    % Reset line properties
    set(ud.hline,'hittest','on','ButtonDownFcn','gnrbcrv(''edit'')')
    % Remove anchor points
    delete(ud.hmark) 
    
case 'endedit'
   disp('endedit')
%    gnrbcrv('editmoveend');
   
   % Get anchor points
   pts = get(ud.hline,'userdata');
   
   if strcmp(mouse(1:3),'alt'),
       pt = get(draw_fig,'CurrentPoint');
       pt = pt(1,1:2);
       pts = insertPoint(pts,pt);      
       set(ud.hmark,'xdata',pts(:,1),'ydata',pts(:,2))
       set(ud.hline,'userdata',pts)
%       return
   end
   
%    set(ud.hline,'UserData',pts,'ButtonDownFcn','gnrbcrv(''edit'')');
   set(draw_fig,'windowbuttonupfcn','gnrbcrv(''editmoveend'')','userdata',ud)
%    set(draw_fig,'userdata',ud)
   set(ud.hline,'ButtonDownFcn','gnrbcrv(''edit'')')
   
otherwise
   error('stats:gnrbcrv:BadCallBack','Invalid call-back.');
end



% Make sure the axis limits don't change
try
    set(ax,'xlim',ud.xlimits,'ylim',ud.ylimits);
catch
end


function [X,Y,Z] = plnlinex(pts,xx,yy,zz)

if isempty(zz), zz=0; end

M = view;
M = M(1:3,1:3);
up = [[1 0 0]; [0 1 0]; [0 0 0]];
M = inv(M)*up;
x = M(1,:) + xx;
y = M(2,:) + yy;
z = M(3,:) + zz;
x = [x pts(1,1) pts(2,1)];
y = [y pts(1,2) pts(2,2)];
z = [z pts(1,3) pts(2,3)];


num = det([ones(1,4); x(1:4); y(1:4); z(1:4)]);
den = det([[ones(1,3); x(1:3); y(1:3); z(1:3)] [0; (x(5)-x(4)); (y(5)-y(4)); (z(5)-z(4))]]);

t = -num/den;

X = x(4) + (x(5)-x(4))*t ;
Y = y(4) + (y(5)-y(4))*t ;
Z = z(4) + (z(5)-z(4))*t ;

% delete(h)

function ind = getClosestPt(h)
x = get(h,'xdata');
y = get(h,'ydata');
z = get(h,'zdata');
if isempty(z), z = x*0; end
n = length(x);
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


function pt = getCurPt(ax)
pt = get(ax,'CurrentPoint');   
pt = pt(1,1:2); 

% Snap to grid if grid on
xgridOn = strcmpi(get(ax,'xgrid'),'on');
ygridOn = strcmpi(get(ax,'ygrid'),'on');
if xgridOn,
    xtick = get(ax,'xtick');
    [m,i] = min(abs(xtick - pt(1)));
    pt(1) = xtick(i);
end
if ygridOn,
    ytick = get(ax,'ytick');
    [m,i] = min(abs(ytick - pt(2)));
    pt(2) = ytick(i);
end


function pts = insertPoint(pts,pt)
if size(pts,2) == 2, 
    pts(:,3) = 0;     
end
pt(:,3) = 0; 
for i=1:size(pt,1),
    d = sqrt( (pts(:,1)-pt(i,1)).^2 + (pts(:,2)-pt(i,2)).^2 + (pts(:,3)-pt(i,3)).^2 );
    if d ~= 0 % Don't add point if it's already on the line
        n = length(d);
        [m,ind] = min(d);
        
        if ind == 1
            % If we're closest to the first point
            % insert just after it
            pts = [pts(1,:); pt(i,:); pts(2:end,:)];
            return
        elseif ind == n
            % If we're closest to the last point
            % insert just before it
            pts = [pts(1:end-1,:); pt(i,:); pts(end,:)];
            return
        end
        
        % If we're closest to an interior point
        A = pts(1:ind,:);
        B = [pts(1:ind-1,:); pt(i,:)];

        dA = diff(A,1,1);
        dB = diff(B,1,1);
        LA = sum(sqrt(dA(:,1).^2 + dA(:,2).^2 + dA(:,3).^2));
        LB = sum(sqrt(dB(:,1).^2 + dB(:,2).^2 + dA(:,3).^2));
        if LA<LB,
            % Insert after closest point
            pts = [pts(1:ind,:); pt(i,:); pts(ind+1:end,:)];
        else
            % Insert before closest point
            pts = [pts(1:ind-1,:); pt(i,:); pts(ind:end,:)];
        end
    end
end

function [x,y,z] = splineInterp(x,y,z)

if nargin < 3,
    z = x*0;
end
x = x(:)';
y = y(:)';
z = z(:)';
if length(x)==1,
    return
end
L = sum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));
npts = L/.01;
pts = [x; y; z];
order = 4;
m = size(pts,2);
n = order + m - 2*(order-1);
while n <= order
    order = (order-1);
    n = order + m - 2*(order-1);
end
% Make a uniform knot sequence
knt = [zeros(1,order-1) linspace(0,1,n) ones(1,order-1)];
crv = nrbmak(pts,knt);
t = linspace(0,1,npts);
p = nrbeval(crv,t);
x = p(1,:)';
y = p(2,:)';
z = p(3,:)';



function [xout,yout,zout] = splineInterp0(x,y,z)

if nargin < 3,
    z = x*0;
end

try
    ctrlpts = [x(:)'; y(:)'; z(:)'];
    sp = cscvn(ctrlpts);
    last = sp.breaks(end);
    npts = length(sp.breaks)*80;
    val = linspace(0,last,npts);
    out = fnval(sp,val);
    xout = out(1,:)';
    yout = out(2,:)';
    zout = out(3,:)';
catch
    % If no spline toolbox we can use a crappier spline interpolation
    disp('Error using cscvn');
    [xout,yout,zout] = splineInterp1(x,y,z);
end






function [xout,yout,zout] = splineInterp1(x,y,z)
% SPLINE_INTERP - Interpolate data with cardinal spline
%   SPLINE_INTERP(P,BC), where P is nx3 matrix of coordinates and BC is a
%   2x3 array containing the tangencies at the beginning and end of the
%   curve

if nargin == 2
    z = zeros(size(x(:)));
elseif isempty(z)
    z = zeros(size(x(:)));
end

P = [x(:) y(:) z(:)];

% Determine number of points
n = size(P,1);

if n <= 2
    xout = x;
    yout = y;
    zout = z;
    return
end

% Determine number of tangencies to account for
m = n - 2;

% Our known vector Reduces to -3*(Pn - Pn+2)
Pa = P(1:m,:);
Pb = P(3:n,:);
PP = -3*(Pa - Pb); %;-[0.5*(1-t)*(1+b)*(1-c)*(Pa-Pb)];

% Create Try diagonal matrix (1,4,1)
TM = (diag(4*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));

% Check for boundary conditions
% if nargin == 2,
    % If no BC specified, impose zero curvature at endpoints
    PP = [6*(P(2,:) - P(1,:)); PP; 6*(P(n,:) - P(n-1,:))];
    TM(1,1:3) = [4 2 0];
    TM(n,n-2:n) = [0 2 4];
% else
%     % Use BC as tangencies at endpoints
%     BC = varargin{1};
%     PP = [BC(1,:); PP; BC(2,:)];
%     TM(1,1:2) = [1 0];
%     TM(n,n-1:n) = [0 1];
% end

% Solve for uknown tangencies 
P_dot = TM\PP; %inv(TM)*PP;     

% Set up matricies for solving
CMx = [P(1:n-1,1) P(2:n,1) P_dot(1:n-1,1) P_dot(2:n,1)];
CMy = [P(1:n-1,2) P(2:n,2) P_dot(1:n-1,2) P_dot(2:n,2)];
CMz = [P(1:n-1,3) P(2:n,3) P_dot(1:n-1,3) P_dot(2:n,3)];

% Calculate interpolated points
num_pts = 79;
res = .1;%(n-1)/(num_pts-1);                % Point frequency
j = 1;

% Sort our data points so that it makes sense for plotting
for u=0:res:1,
    for q = 1:n-1,
        M = [CMx*H(u) CMy*H(u) CMz*H(u)];
        MM(j,:,q) = M(q,:);        
    end
    j = j+1;
end
out = [];
for i=1:n-1,
    out = [out; MM(:,:,i)];
end

% Remove repeated points
out(find(diff(sqrt(sum(out.^2,2)))==0),:) = [];


% Determine length of curve
Length = sum(sqrt(sum(diff(out).^2,2)));
% Curve resolution
% ptsPerInch = 100;
% npts = round(Length*ptsPerInch);
% % Respace points evenly
% out = curvspace(out,npts);
% % Insert Control points back into curve if necessary
% out = insertPoint(out,[x(:) y(:)]);
xout = out(:,1);
yout = out(:,2);
zout = out(:,3);



% Calculate the Hermite Basis Functions
function h = H(u),
h = [(1-3*u^2) + 2*u^3; 3*u^2 - 2*u^3; u - 2*u^2 + u^3; u^3 - u^2;];






















% % Determine number of tangencies to account for
% m = n - 2;
% 
% % Our known vector Reduces to -3*(Pn - Pn+2)
% Pa = P(1:m,:);
% Pb = P(3:n,:);
% PP = -3*(Pa - Pb); %;-[0.5*(1-t)*(1+b)*(1-c)*(Pa-Pb)];
% 
% % Create Try diagonal matrix (1,4,1)
% TM = (diag(4*ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
% 
% % Check for boundary conditions
% if nargin == 2,
%     % If no BC specified, impose zero curvature at endpoints
%     PP = [6*(P(2,:) - P(1,:)); PP; 6*(P(n,:) - P(n-1,:))];
%     TM(1,1:3) = [4 2 0];
%     TM(n,n-2:n) = [0 2 4];
% else
%     % Use BC as tangencies at endpoints
% %     BC = varargin{1};
%     PP = [BC(1,:); PP; BC(2,:)];
%     TM(1,1:2) = [1 0];
%     TM(n,n-1:n) = [0 1];
% end
% 
% % Solve for uknown tangencies 
% P_dot = TM\PP; %inv(TM)*PP;     
% 
% % Set up matricies for solving
% CMx = [P(1:n-1,1) P(2:n,1) P_dot(1:n-1,1) P_dot(2:n,1)];
% CMy = [P(1:n-1,2) P(2:n,2) P_dot(1:n-1,2) P_dot(2:n,2)];
% CMz = [P(1:n-1,3) P(2:n,3) P_dot(1:n-1,3) P_dot(2:n,3)];
% 
% % Calculate interpolated points
% num_pts = 79;
% res = .1;%(n-1)/(num_pts-1);                % Point frequency
% j = 1;
% 
% % Sort our data points so that it makes sense for plotting
% for u=0:res:1,
%     for q = 1:n-1,
%         M = [CMx*H(u) CMy*H(u) CMz*H(u)];
%         MM(j,:,q) = M(q,:);        
%     end
%     j = j+1;
% end
% out = [];
% for i=1:n-1,
%     out = [out; MM(:,:,i)];
% end
% 
% % Remove repeated points
% out(find(diff(sqrt(sum(out.^2,2)))==0),:) = [];

