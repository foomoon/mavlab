function [varargout] = mesh2d(node,cnect,hdata,fix)

% [p,t] = mesh2d(node,cnect,hdata,fix)
% 
% 2D unstructured triangualar mesh generation. 
% 
% A 2D unstructured mesh is generated based on a piecewise-linear boundary
% geometry. An iterative method is implemented to optimise mesh quality. 
% General multiply connected domains can be specified.
%
%   node  = [x1 y1; x2 y2; etc], endpoints of boundaries
%
%   cnect = [n1 n2; n2 n3; etc], boundary connectivity via node numbers
%
%   fix   = [x1 y1; x2; y2; etc], nodes included in the mesh at fixed 
%                                 points
%
% The mesh size can be controlled via a size function specified in the
% structure array hdata. A combination of 3 size function types can be
% specified:
%
%   hdata.type = 1            - Point based size function
%              = 2            - Line based size function  
%              = 3            - Constant size function
%   hdata.slope               - Gradient of size function
%   hdata.h0                  - Reference size at (type 1) the point, 
%                               (type 2) the line or (type 3) everywhere
%   hdata.ref = [x,y]         - the point (type 1)
%             = [x1,y1,x2,y2] - the line endpoints (type 2)
%   hdata.lim                 - distance from ref that the size function is 
%                               applied 
%
% Returns the final coordinates of the nodes p, and their triangulation t
% (with a counter-clockwise node ordering).
%
% Example:
%
%   node = [
%           0   0
%           1   0
%           1   1
%           0   1
%           ];
%        
%   cnect = [
%           1   2
%           2   3
%           3   4
%           4   1
%           ];
%     
%   fix = node; 
% 
%   hdata(1).type  = 1;
%   hdata(1).slope = 0.1;
%   hdata(1).ref   = [0.25,0.75];
%   hdata(1).h0    = 0.01;
%   hdata(1).lim   = 2;
% 
%   hdata(2).type  = 2;
%   hdata(2).slope = -0.04;
%   hdata(2).ref   = [0.25,0.25,0.75,0.75];
%   hdata(2).h0    = 0.05;
%   hdata(2).lim   = 2;
% 
%   [p,t] = mesh2d(node,cnect,hdata,fix);


% Based on the code Distmesh.
%
%   [1] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
%       SIAM Review, Volume 46 (2), pp. 329-345, June 2004
%
% Darren Engwirda - 2005.


% ERROR CHECKING
if nargin~=4, error('Wrong number of inputs');  end
if nargout>3, error('Wrong number of outputs'); end

if (size(node,2)~=2) || (size(cnect,2)~=2)
    error('Input dimensions');
end

if max(cnect(:))>size(node,1)
    error('Incorrect indexing in cnect. Not enough nodes.')
end


% % GET USER TO CHECK SETUP
% figure, hold on, axis equal, grid on, title(['Geometry (red) and size functions (blue). ' ...
%                                              'Limits and centres of size functions shown'])
% 
% % Walls
% patch('faces',cnect,'vertices',node,'facecolor','none','edgecolor','r')

% Size function
for k = 1:length(hdata)
    switch hdata(k).type
        case 1      % Point based
            
            % Centre
            xc = hdata(k).ref(1);
            yc = hdata(k).ref(2);
            % Distance limit
            dh = hdata(k).lim;
            % Plot circle at dh
            x  = dh*cos(0:0.01:2*pi) + xc;
            y  = dh*sin(0:0.01:2*pi) + yc;
            
%             plot(xc,yc,'bo',x,y,'b-.')
            
        case 2      % Line based 
            
            % Endpoints
            x1 = hdata(k).ref(1); y1 = hdata(k).ref(2);
            x2 = hdata(k).ref(3); y2 = hdata(k).ref(4);
            % Distance limit
            dh = hdata(k).lim;
            % Line normals
            nx = -(y2-y1)/sqrt((x2-x1)^2+(y2-y1)^2);
            ny =  (x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2);
            % Outer points
            x3 = x1+dh*nx; y3 = y1+dh*ny;
            x4 = x2+dh*nx; y4 = y2+dh*ny;
            x5 = x1-dh*nx; y5 = y1-dh*ny;
            x6 = x2-dh*nx; y6 = y2-dh*ny;
            
%             plot([x1,x2],[y1,y2],'bo-',[x3,x4],[y3,y4],'b-.',[x5,x6],[y5,y6],'b-.', ...
%                  [x5,x3],[y5,y3],'b-.',[x6,x4],[y6,y4],'b-.')
            
        case 3
        otherwise 
            error('Unknown size function type')   
    end
end

% % User response
% switch lower(questdlg('Ok...???','Geometry and Size Function'))
%     case 'yes'
%         close(gcf)
%     case {'no','cancel'}
%         p = []; t = []; return
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GENERATE MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc, 
tic

% Fixed nodes
nfix = size(fix,1);

% Walls
wall = [node(cnect(:,1),:),node(cnect(:,2),:)];

% Wall endpoints
c = wall(:,1); d = wall(:,2); 
e = wall(:,3); f = wall(:,4); 

% Wall parameters
emc = e-c; 
fmd = f-d;
Lsq = emc.^2+fmd.^2;

% Iteration parameters
qlim    = 0.4;   
Fscale  = 1.2; 
dt      = 0.2; 
maxit   = 200; 
restart = 20;
retri   = true;

% Initial nodes
p = [fix; node];

% Progress bar
w = waitbar(0,'Generating Mesh...');

for iter = 1:maxit

    % Number of nodes
    N = size(p,1);

    % ==========================================================
    %                 Delaunay re-triangulation
    % ==========================================================

    if retri
        
        % New triangles
        t = delaunay2d(p);
        
        % Check for internal triangles
        in = inpoly((p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3,node,cnect);
        
        % Nodes in external triangles
        tout = t(~in,:);
        out  = unique(tout(:));
        % Project them onto the closest wall
        for j = 1:length(out)
            % Current node
            cn = out(j);
            % xy
            x = p(cn,1);
            y = p(cn,2);
            % Parametric ratio
            r = ((x-c).*emc+(y-d).*fmd)./Lsq;
            % Limit to within endpoints
            r = min(max(r,0),1);
            % Normal point on each wall
            xn = c + r.*emc;
            yn = d + r.*fmd;
            % Find closest wall
            [dist,k] = min((xn-x).^2+(yn-y).^2);
            % Move node
            p(cn,1) = xn(k);
            p(cn,2) = yn(k);
        end
        t = t(in,:);
        
        % Describe each bar by a unique pair of nodes
        bars = unique([t(:,[1,2]); t(:,[1,3]); t(:,[2,3])],'rows');
        
        % Bar end nodes
        b1 = bars(:,1); b2 = bars(:,2);
        
        % Size function
        h     = sfun(p,hdata);
        h0    = min(h);
        hbars = 0.5*(h(b1)+h(b2))/h0;
        
    end
    
%     % Animate
%     newplot
%     patch('faces',t,'vertices',p,'facecolor','w','edgecolor','g');
%     axis equal, axis off
%     drawnow

    % ==========================================================
    %             Move nodes based on bar lengths
    % ==========================================================
    
    % Bar forces
    barv = p(b1,:)-p(b2,:);         
    L    = sqrt(sum(barv.^2,2));        
    L0   = hbars*Fscale*(sum(L)/sum(hbars));
    F_L  = min(max(L0-L,0),2*L)./L;

    % Sum bar forces to nodes
    Ftot = zeros(N,2);
    for k = 1:size(bars,1)
        % Nodes
        n1 = b1(k); n2 = b2(k);
        % Forces
        F_x = F_L(k)*barv(k,1);
        F_y = F_L(k)*barv(k,2);
        % Nodal forces
        Ftot(n1,1) = Ftot(n1,1) + F_x;
        Ftot(n1,2) = Ftot(n1,2) + F_y;
        Ftot(n2,1) = Ftot(n2,1) - F_x;
        Ftot(n2,2) = Ftot(n2,2) - F_y;
    end
    if nfix>0   % Don't move fixed
        Ftot(1:nfix,:) = 0;
    end
    
    % Update node positions
    pnew = p + dt*Ftot;
    
    % Find nodes that have moved outside
    % the walls
    in  = inpoly(pnew,node,cnect);
    vec = 1:N;
    out = vec(~in);
    
    % Project them onto the closest wall
    for j = 1:length(out)
        % Current node
        cn = out(j);
        % Parametric ratio
        r = ((pnew(cn,1)-c).*emc+(pnew(cn,2)-d).*fmd)./Lsq;
        % Limit to within endpoints
        r = min(max(r,0),1);
        % Normal point on each wall
        xn = c + r.*emc;
        yn = d + r.*fmd;
        % Find closest
        [dist,k] = min((xn-p(cn,1)).^2+(yn-p(cn,2)).^2);
        % Project node
        pnew(cn,1) = xn(k);
        pnew(cn,2) = yn(k);
    end
    p = pnew;
    
    % ==========================================================
    %                     Break criteria
    % ==========================================================
    
    % Mesh quality
    q    = quality(p,t); 
    minq = min(q);
    
    % Movement within step     
    Lnew = sqrt(sum((p(b1,:)-p(b2,:)).^2,2));
    dL   = Lnew-L;
    i    = Lnew>0;
    move = max(norm(dL(i)./Lnew(i),inf),eps);
  
    % Ratio of actual to desired edge length
    r    = Lnew./(h0*hbars);
    maxr = max(r);
    
    % BREAK CRITERIA
    if (minq>qlim) && (move<0.05) && (maxr<2)
        break
    end
    
    % RE-TRIANGULATION CRITERIA
    if (move>0.25) || (minq<0.1)
        retri = true; 
    else
        retri = false;
    end
    
    
    % ==========================================================
    %                    Density control
    % ==========================================================
    
    if iter<maxit
        % Remove nodes if the edge is less 
        % than half what it should be
        less = r<0.5;
        if sum(less)>0
            % Problem nodes
            prob               = false(N,1);
            prob(bars(less,:)) = true;
            prob(1:nfix)       = true;
            % New nodes
            pnew = [fix; p(~prob,:)];
            % Force retri
            retri = true;
        else
            pnew = p;
        end
        % Add nodes if the edge is more 
        % than twice as long as it should be
        more = r>2;
        if sum(more)>0
            % New nodes
            pnew  = [pnew; 0.5*(p(b1(more),:)+p(b2(more),:))];
            % Force retri
            retri = true;
        end
        % Add nodes at the centroids 
        % of low q triangles
        if iter==restart
            low     = q<qlim;
            add     = (p(t(low,1),:)+p(t(low,2),:)+p(t(low,3),:))/3;
            p       = [pnew; add];
            restart = restart+10;
            % Force retri
            retri = true;
        else
            p = pnew;
        end
    end
    
    % Show progress
    if ~mod(iter,5)
        c1 = min(minq/qlim,1);
        c2 = min(2/maxr,1);
        c3 = min(0.05/move,1);
        waitbar((c1+c2+c3)/3,w);
    end

end
close(w)


% Print stats
Cost = struct('Iterations'  ,iter     ,...
              'Time'        ,toc      ,...
              'Triangles'   ,size(t,1),...
              'Nodes'       ,size(p,1),...
              'Mean_quality',mean(q)  ,...
              'Min_quality' ,min(q));         

% CCW orientation  
[p,t] = fixmesh(p,t);          
          
varargout{1} = p;
varargout{2} = t;
varargout{3} = Cost;

if ~nargout,
    % Final size function
    figure, trimesh(t,p(:,1),p(:,2),sfun(p,hdata)), title('Size Function')

    % Final mesh
    figure, patch('faces',t,'vertices',p,'facecolor','w','edgecolor','g'); hold on, axis equal, axis off

    % Walls
    patch('faces',cnect,'vertices',node,'facecolor','none','edgecolor','r'), title('Mesh')
    clc
    Cost
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = quality(p,t)

% Approximate simplex quality
% From [1].

% xy nodes
p1 = p(t(:,1),:); 
p2 = p(t(:,2),:); 
p3 = p(t(:,3),:);

% Approximate quality
d12 = p2-p1;
d13 = p3-p1;
d23 = p3-p2;
A   = abs(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1));
q   = 2*sqrt(3)*A./(sum(d12.^2,2)+sum(d13.^2,2)+sum(d23.^2,2));

return

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,t] = fixmesh(p,t)

% Remove duplicate nodes & setup CCW node ordering
% From [1].

snap = max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx] = unique(round(p/snap)*snap,'rows');
p = p(ix,:);
t = jx(t);

[pix,ix,jx] = unique(t);
t = reshape(jx,size(t));
p = p(pix,:);


t1 = t(:,1); t2 = t(:,2); t3 = t(:,3);

x12 = p(t2,1)-p(t1,1);
y12 = p(t2,2)-p(t1,2);
x13 = p(t3,1)-p(t1,1);
y13 = p(t3,2)-p(t1,2);

% Negative area?
flip = abs(x12.*y13-y12.*x13)<0;

t(flip,[1,2]) = t(flip,[2,1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = sfun(p,hdata)

% Evaluate size function at the points in p.

x = p(:,1);
y = p(:,2);

num = length(hdata);
h   = inf*ones(size(p,1),num);

for k = 1:num
    switch hdata(k).type
        case 1      % Circular size function
            
            % Centre
            xc = hdata(k).ref(1);
            yc = hdata(k).ref(2);
            % Distance to (xc,yc)
            d = sqrt( (x-xc).^2+(y-yc).^2 );
            % Range
            in = d<hdata(k).lim;
            % Size function
            h(in,k) = hdata(k).h0 + hdata(k).slope*d(in);
            
        case 2      % Line based size function 
            
            % Endpoints
            x1 = hdata(k).ref(1); y1 = hdata(k).ref(2);
            x2 = hdata(k).ref(3); y2 = hdata(k).ref(4);
            % Normal points on line
            r  = ((x-x1)*(x2-x1)+(y-y1)*(y2-y1))/( (x2-x1)^2+(y2-y1)^2 );
            r  = min(max(r,0),1);
            xn = x1 + r*(x2-x1);
            yn = y1 + r*(y2-y1);
            % Normal distance
            d = sqrt( (xn-x).^2+(yn-y).^2 );
            % Range
            in = d<hdata(k).lim;
            % Size function
            h(in,k) = hdata(k).h0 + hdata(k).slope*d(in);
            
        case 3      % Constant size
            
            h(:,k) = hdata(k).h0;
  
    end
end

% Take minimum
h = min(h,[],2);

% Can't be negative
if min(h)<=0
    error('Size function cannot be negative'); 
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in = inpoly(p,node,cnect)

% Determine which points in p are inside the polygonal region defined by
% node/cnect via the crossing number test. Faster than inpolygon.

n = size(p,1); tol = sqrt(eps); 

% Setup walls (wall = [x1 y1 x2 y1; etc])
wall = [node(cnect(:,1),:),node(cnect(:,2),:)];
numw = size(wall,1);
mnyp = sum([wall(:,2); wall(:,4)])/(2*numw);

% Sort walls by y values
[w,i] = sort(min(wall(:,[2 4]),[],2));
[w,j] = sort(max(wall(:,[2 4]),[],2));

% Wall endpoints (not sorted)
Cb = wall(i,1); Db = wall(i,2); Eb = wall(i,3); Fb = wall(i,4);
Ct = wall(j,1); Dt = wall(j,2); Et = wall(j,3); Ft = wall(j,4);

% Wall endpoints (sorted - left, bottom, right, top)
x1b = min(Cb,Eb); y1b = min(Db,Fb); x2b = max(Cb,Eb); y2b = max(Db,Fb); 
x1t = min(Ct,Et); y1t = min(Dt,Ft); x2t = max(Ct,Et); y2t = max(Dt,Ft); 

% Endpoint for test lines
X = 1.1*max(node(:,1));


% Test each point in p

in = false(n,1);
i  = 1;
while i<=n
    
    % Current point
    x = p(i,1); 
    y = p(i,2);
    
    % Initialise
    cn = 0; 
    on = false;
    
    % If the point is in the effective lower half
    if y<=mnyp
        
        % Loop through walls bottom to top
        j = 1;
        while j<=numw
            
            % Partition the test
            if y>=y1b(j)
                if y<=y2b(j)
                    if x>=x1b(j)
                        if x<=x2b(j)
                            % True if on wall
                            c = Cb(j); d = Db(j); 
                            e = Eb(j); f = Fb(j);
                            on = on || (abs((f-y)*(c-x)-(d-y)*(e-x))<tol); 
                            if ~on && (y~=y2b(j))
                                % Check crossing
                                ub = ((e-c)*(d-y)-(f-d)*(c-x))/(-(X-x)*(f-d)); 
                                cn = cn + double((ub>-tol) && (ub<(1+tol)));
                            end
                        end
                    elseif y~=y2b(j)
                        % Has to cross
                        cn = cn+1;
                    end
                end
            else
                % Due to the sorting
                break    
            end 
            
            % Wall counter
            j = j+1;
            
        end
        
    % If the point is in the effective top half    
    else
        
        % Loop through walls top to bottom
        j = numw;
        while j>0
            % Partition the test
            if y<=y2t(j)
                if y>=y1t(j)
                    if x>=x1t(j)
                        if x<=x2t(j)
                            % True if on wall
                            c = Ct(j); d = Dt(j); 
                            e = Et(j); f = Ft(j);
                            on = on || (abs((f-y)*(c-x)-(d-y)*(e-x))<tol); 
                            if ~on && (y~=y2t(j))
                                % Check crossing
                                ub = ((e-c)*(d-y)-(f-d)*(c-x))/(-(X-x)*(f-d)); 
                                cn = cn + double((ub>-tol) && (ub<(1+tol)));
                            end
                        end
                    elseif y~=y2t(j)
                        % Has to cross
                        cn = cn+1;
                    end
                end
            else 
                % Due to the sorting
                break    
            end 
            
            % Wall counter
            j = j-1;
            
        end 
        
    end
    
    % Point is inside if cn is odd or if it is on a wall
    in(i) = mod(cn,2) || on;
    
    % Counter
    i = i+1;
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = delaunay2d(x)

% My version of the MATLAB function delaunayn for 2D tesselations
% I modified the section to remove simplicies with ~zero area to 
% improve run time.

[x,idx,jdx] = unique(x,'rows');
[m,n]       = size(x);

if m<n+1,
  error('Not enough unique points to do tessellation.');
end
if any(isinf(x(:)) | isnan(x(:)))
  error('Data containing Inf or NaN cannot be tessellated.');
end
if m==n+1
  t = 1:n+1;
  return;
end

% qhull needs the sum of squares at the end of the points
t = qhullmx([x,sum(x.*x,2)]', 'd ', 'QJ ', 'Pp');


% Deal with ~zero area triangles

t1 = t(:,1); t2 = t(:,2); t3 = t(:,3);

x12 = x(t2,1)-x(t1,1);
y12 = x(t2,2)-x(t1,2);
x13 = x(t3,1)-x(t1,1);
y13 = x(t3,2)-x(t1,2);

% Tolerance
seps = sqrt(eps)*norm(x,'inf');

t = t(abs(x12.*y13-y12.*x13)>seps,:);
t = idx(t);

return
