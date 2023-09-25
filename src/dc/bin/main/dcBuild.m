function [X,Y,Z] = dcBuild(wing,varargin),
% X - OUT RIGHT WING
% Y - UPSTREAM
% Z - UP

global dcDebug


if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  Wing Regenerated\n')
    fprintf(1,'  Span:     %0.3f \n',wing.span)
    fprintf(1,'  Chord:    %0.3f \n',wing.chord)
    fprintf(1,'  Camber:   %0.3f \n',wing.camber)
    fprintf(1,'  Tip:      %0.2f \n',wing.tip)
    fprintf(1,'  Twist:    %0.1f \n',wing.twist)
    fprintf(1,'  Sweep:    %0.1f \n',wing.sweep)
    fprintf(1,'  Dihedral: %0.1f \n',wing.dihedral)
    fprintf(1,'  Edge Ref: %0.0f \n',wing.edgeref)    
end

if isfield(wing,'ledge')
    % Then we're using and old MAV file
    [X,Y,Z] = build_wing3(wing);    
    if dcDebug              
       fprintf(1,'\n  Using Old Format\n')  
    end
    return
end

% Get Planform Leading Edge
LedgeX = wing.planform.Ledge.x(:)';
LedgeY = wing.planform.Ledge.y(:)';
try
    LedgeZ = wing.planform.Ledge.z(:)';
catch
    LedgeZ = zeros(size(LedgeX));
end

% Get Planform Trailing Edge
TedgeX = wing.planform.Tedge.x(:)';
TedgeY = wing.planform.Tedge.y(:)';
try
    TedgeZ = wing.planform.Tedge.z(:)';
catch
    TedgeZ = zeros(size(LedgeX));
end

% Shift wing if necessary ( Trailing Edge root should be at [0,0,0] )
xshift = -LedgeX(1);
yshift = -TedgeY(1);
zshift = -LedgeZ(1);
LedgeX = LedgeX + xshift; TedgeX = TedgeX + xshift;
LedgeY = LedgeY + yshift; TedgeY = TedgeY + yshift;
LedgeZ = LedgeZ + zshift; TedgeZ = TedgeZ + zshift;

% Spline Leading/Trailing Edge
if nargin > 1
    xnpts = varargin{1};
else
    xnpts = wing.span*10/2;
end

if ~isempty(xnpts)
%     xii = linspace(min(LedgeX),max(LedgeX),xnpts);
%     xii = cosdist(xnpts)*max(LedgeX);
    theta = linspace(pi/2,pi,xnpts);
    xii = -max(LedgeX)*cos(theta);
    xii(1) = 0;

    % Interpolate leading/trailing edges onto grid
    method = 'cubic';
    LedgeY = interp1(LedgeX,LedgeY,xii,method); % Leading  Edge
    TedgeY = interp1(TedgeX,TedgeY,xii,method); % Trailing Edge
    LedgeZ = interp1(LedgeX,LedgeZ,xii,method); % Leading  Edge
    TedgeZ = interp1(TedgeX,TedgeZ,xii,method); % Trailing Edge
    LedgeX = xii;
    TedgeX = xii;
end

% Find wing tip
Xmax  = max([max(abs(LedgeX)) max(abs(TedgeX))]);
Zmax = LedgeZ(end); if Zmax == 0, Zmax = 1; end

% Find root chord size
Yroot = LedgeY(1); 

% Normalize planform data
LedgeX = LedgeX/Xmax; %TedgeX = TedgeX/Xmax;
LedgeY = LedgeY/Yroot; TedgeY = TedgeY/Yroot;
LedgeZ = LedgeZ/Zmax; %TedgeZ = TedgeZ/Zmax;

% % DEBUG
% figure; 
% line(LedgeX,LedgeY,LedgeZ,'color','r');
% line(TedgeX,TedgeY,TedgeZ,'color','b');
% axis equal

% Chord Widths
Cwidth = LedgeY - TedgeY;


% FORMAT AIRFOIL DATA
[AirfoilX,AirfoilY] = buildAirfoils(wing,LedgeX,LedgeY,Cwidth);


% CALCULATE MAX CAMBER REFERENCE DISTANCE
n = size(AirfoilX,1);
maxZ = max(max(AirfoilY));
maxPerFoil = max(AirfoilY,[],1);
maxOffset = (maxZ-maxPerFoil);
maxOffset = maxOffset(ones(n,1),:);

% SCALE FACTORS
scaleX = wing.span/2;
scaleY = wing.chord;
scaleZ = wing.chord;

% ROTATIONS
phi =   -(pi/180)*wing.twist;
gamma = -(pi/180)*wing.sweep;
theta =  (pi/180)*wing.dihedral;

% TWIST
[AirfoilX,AirfoilY] = twist(LedgeX,Cwidth,AirfoilX,AirfoilY,phi);

% SWEEP
LedgeY = LedgeY + LedgeX*tan(gamma);

% DIHEDRAL
LedgeZ = LedgeZ*wing.tip + LedgeX*tan(theta)*scaleX;

% ASSEMBLE
X = LedgeX(ones(n,1),:) + 0;
Y = LedgeY(ones(n,1),:) - AirfoilX; 
Z = LedgeZ(ones(n,1),:) + (AirfoilY + maxOffset*(1-wing.edgeref))*scaleZ;
% X = X - max(max(X));

% APPLY SCALE FACTOR (FINAL PART OF VARIABLE CHANGE HACK)
X = X*scaleX;
Y = Y*scaleY;

% MIRROR SURFACE ABOUT X = 0
if (isfield(wing,'mirror') & wing.mirror) | ~isfield(wing,'mirror')
    X = [-fliplr(X) X(:,2:end)];
    Y = [fliplr(Y) Y(:,2:end)];
    Z = [fliplr(Z) Z(:,2:end)];
end

% NEED TO TRANSPOSE DATA
X = X';
Y = Y';
Z = Z';



if nargout <= 1,
    wing.X = X;
    wing.Y = Y;
    wing.Z = Z;
    X = wing;
end





function [AirfoilX,AirfoilY] = buildAirfoils(wing,LedgeX,LedgeY,Cwidth)
AirfoilX = wing.foils.x;
AirfoilY = wing.foils.y;

n = size(AirfoilX,1);
m = length(Cwidth);

Ny = linspace(0,1,n);
LedgeX(LedgeX<0) = 0;
[xi,yi] = meshgrid(LedgeX,Ny);

% DETERMINE TYPE OF AIRFOIL BLENDING FROM INDEX
buildType = length(wing.ind);
buildType = min(buildType,2);

% BLEND AIRFOILS
switch buildType,
    case 0,   % This method is similar to a wing stamped out of a curved plate
        LedgeY = LedgeY - LedgeY(1);
        x0 = -LedgeY;
        x1 = -(LedgeY - Cwidth);
        [AirfoilX,AirfoilY] = rolledPlate(AirfoilX(:,1),AirfoilY(:,1),x0,x1,n);
        Xoffset = min(AirfoilX,[],1);
        Yoffset = AirfoilY(1,:); 
        AirfoilX = AirfoilX - Xoffset(ones(size(AirfoilX,1),1),:);
        AirfoilY = AirfoilY - Yoffset(ones(size(AirfoilX,1),1),:);
    case 1,   % Single Airfoil Case
        [AirfoilX,AirfoilY] = blend([0 1],AirfoilX(:,[1 1]),AirfoilY(:,[1 1]),xi,yi);
        AirfoilX = AirfoilX.*Cwidth(ones(n,1),:);
        AirfoilY = AirfoilY.*Cwidth(ones(n,1),:);
    case 2,   % Multi Airfoils
        [AirfoilX,AirfoilY] = blend(wing.ind,AirfoilX,AirfoilY,xi,yi);
        AirfoilX = AirfoilX.*Cwidth(ones(n,1),:);
        AirfoilY = AirfoilY.*Cwidth(ones(n,1),:);
end

maxCamber = max(AirfoilY(:,1));
maxCamber = max(eps,maxCamber);
AirfoilY = AirfoilY/maxCamber*wing.camber; 



function T = Trotate(theta),
% ROTATION MATRIX
T = [cos(theta) sin(theta); -sin(theta) cos(theta)];


function [X,Y] = twist(LX,LY,AirfoilX,AirfoilY,phi),
% SPANWISE TWIST (THETA)
theta = LX*phi;

% SHIFT TO ROTATE ABOUT QUARTER CHORD
Cq = (1/4)*LY(ones(1,size(AirfoilX,1)),:);
AirfoilX = AirfoilX - Cq;

% BUILD ROTATION MATRIX
D0          = zeros(1,2*length(theta));
D0(1:2:end) = cos(theta);
D0(2:2:end) = cos(theta);
D1          = zeros(1,length(D0)-1);
D1(1:2:end) = sin(theta);
D2          = zeros(1,length(D0)-1);
D2(1:2:end) = -sin(theta);
T           = diag(D0) + diag(D1,1) + diag(D2,-1);

% ORGANIZE AIRFOIL DATA ([X1; Y1]; [X2; Z2]...)
U             = zeros(2*size(AirfoilX,2),size(AirfoilX,1));
U(1:2:end,:)  = AirfoilX';
U(2:2:end,:)  = AirfoilY';

% SOLVE
F = T*U;

% SHIFT BACK
X = F(1:2:end,:)' + Cq;
Y = F(2:2:end,:)';


function [Xnew,Ynew] = rolledPlate(X,Y,x0,x1,n),

for i=1:length(x0),
    Xnew(:,i) = linspace(x0(i),x1(i),n)';
    Ynew(:,i) = spline(X,Y,Xnew(:,i));
end


function [yi,zi] = blend(ind,y,z,xi,yi)

sy = size(y);
sz = size(z);

if length(ind) ~= sy(2)
    error('The length of first arguement must be the same as the number of columns in the 2nd and 3rd');
end
if ~~sum(sy - sz) 
    error('X and Y must be the same size');
end
if length(ind) < 2
    error('At least two points are required to blend');
end

x = repmat(sort(ind(:))',sy(1),1);

extrapval = 0; % Value used where data is outside of interpolation region
if sy(2) > 2
    zi = interp2(x,y,z,xi,yi,'cubic',extrapval);
else
    zi = interp2(x,y,z,xi,yi,'linear',extrapval);
end

