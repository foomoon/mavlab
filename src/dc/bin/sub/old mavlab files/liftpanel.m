function [cl,cm,cd,LD] = liftpanel(z,alfa,Re,varargin),
% LIFTPANEL 2D Panel Method Solver
%   LIFTPANEL is a 2D viscous panel method solver.
%
%   [cl,cm,cd,LD]=liftpanel(airfoil,alfa,Re)
%   
%   Input
%       airfoil:  (2D-array) airfoil coordinates
%       alpha:    (scaler) angle of attack (degrees)
%       Re:       (scaler) Reynold's number
%
%   Output
%       cl:       (scaler) coefficient of lift
%       cm:       (scaler) coefficient of pitching moment
%       cd:       (scaler) coefficient of drag
%       LD:       (scaler) Lift/Drag
%
% See also: MAVLAB


if nargin == 3,
    plotON = 0;
elseif nargin > 3,
    plotON = varargin{1};
else
    plotON = 0;
end

nbp=length(z)-1;

%turn viscous solver on (1) or off (0)
viscous=1;

%call solver
[cl,cm,cd,LD]=solver(z,alfa,Re,nbp,viscous,plotON);

% hold off;

return


%%%% FUNCTION %%%%

function [cl,cm,cd,LD]=solver(z,alfa,Re,nbp,viscous,varargin)

%SOLVER
%Essentially keeps the books, calls vortex and sets up the viscous portion
%Modified from Pablo

%inviscid portion

clcm = vortex(z,alfa,varargin{1});

cl = clcm(1);
cm = clcm(2);
ue = clcm(3:nbp+3);
xcp =clcm(nbp+4);
cpmin = clcm(nbp+5);
cpmax = clcm(nbp+6);
 
cd=0;

%viscous portion
if viscous==1
    
    plotbl=0;
    
    % Determine the stagnation point
    isp = stagnation_point(ue);
    ue(isp) = 0;
    
    % Boundary layer discretization
    nup = isp; nlo = nbp+2-isp;
    
    chordL = 1;
    zup = z(isp:-1:1,:)/chordL;       % Normalized length
    zlo = z(isp:nbp+1,:)/chordL;
    
    Vzero = 1;
    ueup = -ue(isp:-1:1)/Vzero;       % Normalized velocity
    uelo = ue(isp:nbp+1)/Vzero;
    
    resup = solvebl(Re,zup,nup,ueup,plotbl,1);
    reslo = solvebl(Re,zlo,nlo,uelo,plotbl,2);
    cd = sy(resup(1),resup(2),resup(3),reslo(1),reslo(2),reslo(3));
end

if viscous==1
    LD=cl./cd;
else
    LD=0;
end

return

%%%% FUNCTION %%%%

function d = dist(z,i,j);

%distance between two points

d = sqrt((z(i,1)-z(j,1)).^2+(z(i,2)-z(j,2)).^2);

return

%%%% FUNCTION %%%%

function res = vortex(za,alfa,varargin);

% This program finds and plots the pressure distribution
% on an airfoil by representing the surface as a finite
% number of linear strength vortex panels.
% (Neumann Boundary condition V.n = 0)
% Airfoils are taken from the Naca 4 digits library.

% Input data
% za is an array containing the airfoil panels coordinates
% alfa is the angle of attack expressed in degrees

if nargin > 2,
    plotON=varargin{1}; %turns plot on/off
else
    plotON = 0;
end

nbp = max(size(za))-1; % number of panels
chord = 1;
Vzero = 1;
alfar = pi.*alfa./180;

% Turn z into an array of complex number
z = za(:,1)+i*za(:,2);

% Change z to clockwise
z = z(nbp+1:-1:1);

% Collocation points
m =(z(1:nbp)+z(2:nbp+1))/2;

% Panel angle
th = imag(log(z(2:nbp+1)-z(1:nbp)));

% Free stream normal velocity component

RHSi(1:nbp,1) = cos(alfar)*sin(th(1:nbp))-sin(alfar)*cos(th(1:nbp));

% Influence matrix

% convert collocation pt to panel cs
xzt = m*ones(1,nbp)-ones(nbp,1)*z(1:nbp).';
xt = real(xzt);
zt = imag(xzt);

xz2t = diff(z);
x2t = real(xz2t);
z2t = imag(xz2t);

cth = ones(nbp,1)*cos(th).';
sth = ones(nbp,1)*sin(th).';

X = xt.*cth+zt.*sth;
Z = -xt.*sth+zt.*cth;
X2 = x2t.*cos(th)+z2t.*sin(th);

% compute r1,r2 and th2-th1

mii = m*ones(1,nbp);
zjj = ones(nbp,1)*z(1:nbp).';
zjjp1 = ones(nbp,1)*z(2:nbp+1).';

r1 = abs(zjj-mii);
r2 = abs(zjjp1-mii);

angle = imag(log((zjjp1-mii)./(zjj-mii))) ;
angle = mod(angle-pi,2*pi)-pi;
tmp = X2(:);
X2mat = ones(nbp,1)*tmp';

th2mth1 = angle./(2*pi*X2mat);
RR = log(r2./r1)./(2*pi*X2mat);
u2l = (Z.*RR+X.*th2mth1);
u1l = -(u2l-X2mat.*th2mth1);
cnst = 1/(2*pi);
TMP = 1/(2*pi)-Z.*th2mth1;
w1l = -TMP+(X2mat-X).*RR;
w2l = TMP+X.*RR;

tmp = diag(u1l) + 0.5*( diag(X)-X2  )./X2;
u1l = u1l - diag( tmp );
tmp = diag(u2l) - 0.5*diag(X)./X2;
u2l = u2l - diag(tmp);
tmp = diag(w1l);
w1l = w1l - diag(tmp + 1/(2*pi));
tmp = diag(w2l);
w2l = w2l - diag(tmp - 1/(2*pi));

% Velocity in global cs

ca = ones(nbp,1)*cos(-th)';
sa = ones(nbp,1)*sin(-th)';

u1 = u1l.*ca+w1l.*sa;
u2 = u2l.*ca+w2l.*sa;
w1 = -u1l.*sa+w1l.*ca;
w2 = -u2l.*sa+w2l.*ca;

% Influence matrix coefficient

CA = cos(th)*ones(1,nbp);
SA = sin(th)*ones(1,nbp);
Aij = zeros(nbp,nbp+1);
Bij = Aij;
Aij(:,1:nbp) = -u1.*SA + w1.*CA;
Bij(:,1:nbp) =  u1.*CA + w1.*SA;
Aij(:,2:nbp+1) = Aij(:,2:nbp+1) - u2.*SA + w2.*CA;
Bij(:,2:nbp+1) = Bij(:,2:nbp+1) + u2.*CA + w2.*SA;

% Add a wake panel with a constant-strength vortex

% Infinite wake point 

d1n = dist(za,nbp,nbp+1);
d12 = dist(za,2,1);
xP = (za(nbp,1)-za(nbp+1,1))./d1n + (za(2,1)-za(1,1))./d12;
yP = (za(nbp,2)-za(nbp+1,2))./d1n + (za(2,2)-za(1,2))./d12;
dPo = sqrt(xP.^2+yP.^2);
zi1 = ((za(1,1)+za(nbp+1,1))/2-chord.*2000.*xP./dPo);
zi2 = ((za(1,2)+za(nbp+1,2))/2-chord.*2000.*yP./dPo);

zi = zi1 + i*zi2;
zte = z(1);

d1 = zte*ones(nbp,1) - m;
d2 = zi*ones(nbp,1) - m;

angle = imag(log(d1./d2)) ;
angle = mod(angle-pi,2*pi)-pi;

r1or2 = abs(d1)./abs(d2);

u = 1/(2*pi)*angle;
w = -1/(2*pi)*log(r1or2);

% transfer to global cs

ca = cos(-th);
sa = sin(-th);

ug = u.*ca+w.*sa;
wg = -u.*sa+w.*ca;

% find the tangential component : 
CA = cos(th);
SA = sin(th);

Aw = -ug.*SA + wg.*CA;

Aij = [Aij Aw];

% Kutta condition

Aij(nbp+1,1)     = 1;
Aij(nbp+2,nbp+1) = 1;

RHSi(nbp+1,1) = 0;
RHSi(nbp+2,1) = 0;

% Solve

gamma = Aij\RHSi;

gamma = gamma(1:nbp+1,1);

% Compute velocity

vel =  Bij*gamma;
velocity = vel + cos(alfar)*cos(th)+sin(alfar)*sin(th);

% turn z back to anti-clockwise
z = z(nbp+1:-1:1); 

% turn velocity to anti-clockwise
velocity = -velocity(nbp:-1:1);

% Velocity at nodes
Qtj = (velocity(1:nbp-1) + velocity(2:nbp))./2;

dz12 = abs(z(1)-z(2));
dz23 = abs(z(3)-z(2));
dznbp1 = abs(z(nbp)-z(1));
dznbpnbpm1 = abs(z(nbp)-z(nbp-1));

Qtj1 = Qtj(1)+dz12.*(Qtj(1)-Qtj(2))./dz23;
Qtj2 = Qtj(nbp-1)+dznbp1.*(Qtj(nbp-1)-Qtj(nbp-2))./dznbpnbpm1;

v(1) = (Qtj1-Qtj2)./2;
v(2:nbp) = Qtj;
v(nbp+1)= - v(1);

% compute cp
cp = 1-velocity.^2/Vzero.^2;

% Graphics

cpmax = max(cp); cpmin = min(cp);
%if cponly==0
  zplot = real(z)-i*(cpmax-cpmin)*imag(z);
  
  if plotON,
%     plot(zplot,'k');hold on; 
    axis('ij'); 
    axis([-0.1 1.1 cpmin cpmax]); 
  end
   
  %end
if plotON,
    plot(real(m),cp,'g'); hold on;
end

% lift coefficient

Fx = 0;
Fy = 0;
cm = 0;
cmle = 0;

for jj=1:nbp
  fxj = -cp(jj)*imag(z(jj+1)-z(jj));
  fyj = cp(jj)*real(z(jj+1)-z(jj));
  Fx = Fx + fxj;
  Fy = Fy + fyj;
  cm = cm + fxj*imag(m(jj)) - fyj*(real(m(jj))-chord/4);
  cmle = cmle + fxj*imag(m(jj)) - fyj*real(m(jj));
end

cl = Fy*cos(alfar) - Fx*sin(alfar);

xcp = -cmle/Fy;

if abs(cl)<0.001
    cl = 0;
    cm = 0;
    xcp = 0;
else
    cl = floor(10000*cl)/10000;
    cm = floor(10000*cm)/10000;
end;

xcp = floor(100*xcp)/100;

res(1) = cl;
res(2) = cm;
res(3:nbp+3) = v;
res(nbp+4) = xcp;
res(nbp+5) = cpmin;
res(nbp+6) = cpmax;

%end

%THE FOLLOWING FUNCTIONS ARE REQUIRED FOR THE VISCOUS PORTION ONLY
%AND CAN BE REMOVED FROM THE INVISCID CODE

%%%% FUNCTION %%%%

function isp = stagnation_point(Qtj);

isp = 2; 

if Qtj(isp)>0
   while Qtj(isp)>0  
       isp = isp+1;
   end;
elseif Qtj(isp)<0
   while Qtj(isp)<0  
       isp = isp+1;
   end;
end;

ispu = isp -1;
ispl = isp;

if abs(Qtj(ispu)) < abs(Qtj(ispl)) 
   isp = ispu;
else
   isp = ispl;
end; 

return

%%%% FUNCTION %%%%

function res = solvebl(Re,z,n,ue,plotcp,side);

nbp2 = 200;  % bl discretization

% Arc length
ss(1,1) = 0;
for ii=2:n
  ss(ii,1) = ss(ii-1,1)+dist(z,ii,ii-1);
end;

sTE = ss(n,1);

% compute the boundary layer up to xe

Coff = 0.98;

nm = floor(0.5*n);
se = spline(z(nm:n,1),ss(nm:n),Coff); 

s = 0:se/nbp2:se;

% Velocity

spues = spline(ss,ue);
ues = ppval(spues,s);

% detect if there is some ue < 0
% which means a prb with the spline interpolation

in = find(ues <0);

if ~isempty(in)

  % Add some more points near the LE
  ue2(1) = ue(1);
  ue2(2) = 0.5*(ue(2)+ue(1));
  ue2(3) = ue(2);
  ue2(4) = 0.5*(ue(3)+ue(2));
  ue2(5) = ue(3);
  ue2(6) = 0.5*(ue(4)+ue(3));
  ue2(7:n+3) = ue(4:n);

  ss2(1) = ss(1);
  ss2(2) = 0.5*(ss(2)+ss(1));
  ss2(3) = ss(2);
  ss2(4) = 0.5*(ss(3)+ss(2));
  ss2(5) = ss(3);
  ss2(6) = 0.5*(ss(4)+ss(3));
  ss2(7:n+3) = ss(4:n);

  % re-spline
  spues = spline(ss2,ue2);
  ues = ppval(spues,s);

end;

% x coordinate

spx = spline(ss,z(:,1));

%plot(s,ues,'k');
%hold on;
%plot(ss,ue,'o');
%pause

ue = ues;

n = nbp2+1;

% x coordinate

spx = spline(ss,z(:,1));

% velocity gradient at nodes

v1= ue(1);  v2 = ue(2);  v3 = ue(3);
x1= s(1); x2 = s(2); x3 = s(3);
gamma = 1./(x3-x2)*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(1) = gamma*(x1-x2) + fac;

if dueds(1) < 0
  dueds(1) = (v2-v1)/(x2-x1);
end;

v1=ue(1:n-2); v2=ue(2:n-1); v3=ue(3:n);
x1=s(1:n-2);x2=s(2:n-1);x3=s(3:n);
gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(2:n-1) = gamma.*(x2-x1) + fac;

v1= ue(n-2);  v2 = ue(n-1);  v3 = ue(n);
x1= s(n-2); x2 = s(n-1); x3 = s(n);
gamma = 1./(x3-x2)*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
fac = (v2-v1)./(x2-x1);
dueds(n) = gamma*(2*x3-x1-x2) + fac;

%--------Laminar boundary layer

lsep = 0; trans=0; endofsurf=0;

theta(1) = sqrt(0.075/(Re*dueds(1)));
i = 1;

while lsep ==0 & trans ==0 & endofsurf ==0

  lambda = theta(i).^2*dueds(i)*Re;

  % test for laminar separation
  if lambda < -0.09 
    lsep = 1;
    itrans = i;
    break; 
  end;

  H(i) = fH(lambda);
  L = fL(lambda);

  cf(i) = 2*L./(Re*theta(i));
  if i>1, cf(i) = cf(i)./ue(i); end;
  i = i+1;

  % test for end of surface
  if i> n endofsurf = 1; itrans = n; break; end;  

  K = 0.45/Re;
  xm = (s(i)+s(i-1))/2;
  dx = (s(i)-s(i-1));
  coeff = sqrt(3/5);

  f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
  f2 = ppval(spues,xm);            f2 = f2^5;
  f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

  dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
  theta(i) = sqrt((theta(i-1).^2*ue(i-1).^6 + dth2ue6)./ue(i).^6);

  % test for transition
  rex = Re*s(i)*ue(i);
  ret = Re*theta(i)*ue(i);
  retmax = 1.174*(rex^0.46+22400*rex^(-0.54));
  if ret>retmax 
    trans = 1; 
    itrans = i;
  end;

end;

%-------- Transition 

transorlamsep = 0;
transloc = 1;
tsep = 0;

if itrans < n

  if trans == 1

     uei = ue(i); thi = theta(i); si = s(i); duedsi = dueds(i);
     ueim1 = ue(i-1); thim1 = theta(i-1); sim1 = s(i-1); duedsim1 = dueds(i-1);

     % Find f(x) at i and i-1

     fxi = ret - retmax;  % already computed 
 
     rex = Re*sim1*ueim1;
     ret = Re*thim1*ueim1;
     retmax = 1.174*(rex^0.46+22400*rex^(-0.54));  
     fxim1 = ret - retmax;

     % Fit a linear function and find the root

     st = sim1 - fxim1/((fxi-fxim1)/(si-sim1));

     transorlamsep = 1;
     transloc = 100*ppval(spx,st);

     % Find the value of theta, and H at st using thwaites

     uet = ppval(spues,st);

     v1=ueim1; v2=uet; v3=uei;
     x1=sim1;  x2=st;  x3=si;
     gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
     fac = (v2-v1)./(x2-x1);
     duedst = gamma.*(x2-x1) + fac;

     xm = (st+sim1)/2;
     dx = (st-sim1);

     f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
     f2 = ppval(spues,xm);            f2 = f2^5;
     f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

     dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
     thetat  = sqrt((thim1.^2*ueim1.^6 + dth2ue6)./uei.^6);

     lambdat = thetat.^2*duedst*Re;
     Ht = fH(lambdat);

     if Ht < 1.1
        Ht = 1.2;
     end;

     if Ht > 2  % to avoid turbulent separation just after transition
        Ht = 2;
     end;

     % Find the value of theta, and H at i using head

     y(1) = thetat;
     y(2) = H1ofH(Ht);
     dx = s(i) - st;

     y = runge(dx,y,Re,uet,duedst,ue(i),dueds(i));

     theta(i) = y(1);
     H(i) = HofH1(y(2));
     rtheta = Re*ue(i)*theta(i);
     cf(i) = cfturb(rtheta,H(i));     

elseif lsep == 1

     uei = ue(i); thi = theta(i); si = s(i); duedsi = dueds(i);
     ueim1 = ue(i-1); thim1 = theta(i-1); sim1 = s(i-1); duedsim1 = dueds(i-1);

     % Find f(x) at i and i-1

     fxi = thi.^2*duedsi*Re+0.09;
     fxim1 = thim1.^2*duedsim1*Re+0.09;
 
     % fit a linear function and find the root

     st = sim1 - fxim1/((fxi-fxim1)/(si-sim1));
 
     transorlamsep = 2;
     transloc = 100*ppval(spx,st);

     % Find the value of theta, and H at st using thwaites

     uet = ppval(spues,st);

     v1=ueim1; v2=uet; v3=uei;
     x1=sim1;  x2=st;  x3=si;
     gamma = 1./(x3-x2).*((v3-v1)./(x3-x1)-(v2-v1)./(x2-x1));
     fac = (v2-v1)./(x2-x1);
     duedst = gamma.*(x2-x1) + fac;

     xm = (st+sim1)/2;
     dx = (st-sim1);

     f1 = ppval(spues,xm-coeff*dx/2); f1 = f1^5;
     f2 = ppval(spues,xm);            f2 = f2^5;
     f3 = ppval(spues,xm+coeff*dx/2); f3 = f3^5;

     dth2ue6 = K*dx/18*(5*f1+8*f2+5*f3);
     thetat  = sqrt((thim1.^2*ueim1.^6 + dth2ue6)./uei.^6);

     lambdat = thetat.^2*duedst*Re;
     Ht = fH(lambdat);

     if Ht < 1.1
        Ht = 1.2;
     end;

     if Ht > 2  % to avoid turbulent separation just after transition
        Ht = 2;
     end;

     % Find the value of theta, and H at i using head

     y(1) = thetat;
     y(2) = H1ofH(Ht);
     dx = s(i) - st;

     y = runge(dx,y,Re,uet,duedst,ue(i),dueds(i));

     theta(i) = y(1);
     H(i) = HofH1(y(2));
     rtheta = Re*ue(i)*theta(i);
     cf(i) = cfturb(rtheta,H(i));     
     
end;

%--------TURBULENT BL

  tsep = 0;
  
  i = i+1;
  
  while endofsurf == 0 & tsep ==0;

    y = runge(s(i)-s(i-1),y,Re,ue(i-1),dueds(i-1),ue(i),dueds(i));
 
    theta(i) = y(1);
    H(i) = HofH1(y(2));

    if H(i) == 3 % which is actually a flag 
        tsep = 100*ppval(spx,s(i));
        i = i-1;
    end;

    rtheta = Re*ue(i)*theta(i);
    cf(i) = cfturb(rtheta,H(i));

    i = i+1;

    if i>n endofsurf = 1; end;

  end;

end;

% plotcp=0;
if plotcp == 1

    deltas = H(1:i-1).*theta(1:i-1);

    if side ==1
  
      figure;plot(s(1:i-1),theta(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Momentum Thickness Theta');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),deltas);grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Displacement Thickness delta star');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),H(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Shape Factor H');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),cf(1:i-1));grid;
      h = title('Upper Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Skin Friction Coefficient Cf');set(h,'Fontsize',[14]);

    elseif side == 2

      figure;plot(s(1:i-1),theta(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Momentum Thickness Theta');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),deltas);grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Displacement Thickness delta star');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),H(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Shape Factor H');set(h,'Fontsize',[14]);

      figure;plot(s(1:i-1),cf(1:i-1));grid;
      h = title('Lower Side');set(h,'Fontsize',[14]); 
      h = xlabel('Arc Length s');set(h,'Fontsize',[14]); 
      h = ylabel('Skin Friction Coefficient Cf');set(h,'Fontsize',[14]);

    end;

end;

res(1) = theta(i-1);
res(2) = H(i-1);
res(3) = ue(i-1);
res(4) = transorlamsep;  % =0 is fully laminar, =1 if transition, =2 if laminar separation
res(5) = floor(100*transloc)/100;
res(6) = floor(100*tsep)/100; % if =0, means that there is no turbulent separation

return

%%%% FUNCTION %%%%

function H = fH(lambda);

if lambda < 0

  if lambda==-0.14
    lambda=-0.139;
    disp('H(lambda) : Lambda = -0.14 -> Lambda = -0.139'); 
  end;

  H = 2.088 + 0.0731./(lambda+0.14);

elseif lambda >= 0 

  H = 2.61 - 3.75*lambda + 5.24*lambda.^2;

end;

return

%%%% FUNCTION %%%%

function L = fL(lambda);

if lambda < 0

  if lambda==-0.107
    lambda=-0.106;
    disp('l(lambda) : Lambda = -0.107 -> Lambda = -0.106');
  end;

  L = 0.22 + 1.402*lambda + (0.018*lambda)./(lambda+0.107);

elseif lambda >= 0 

  L = 0.22 + 1.57*lambda - 1.8*lambda.^2;

end;

return

%%%% FUNCTION %%%%

function y=H1ofH(H);

if H <1.1
   disp('H < 1.1 !  -> H1 = 16');
   y = 16;
else   
  if H <= 1.6       
    y = 3.3 + 0.8234*(H-1.1).^(-1.287);  
  else
    y = 3.3 + 1.5501*(H-0.6778).^(-3.064);
  end;
end;

return

%%%% FUNCTION %%%%

function H=HofH1(H1);

if H1 <= 3.32
  H = 3;
elseif H1 < 5.3
  H = 0.6778 + 1.1536*(H1-3.3).^(-0.326);
else
  H = 1.1 + 0.86*(H1-3.3).^(-0.777);
end

return

%%%% FUNCTION %%%%

function ynp1 = runge(dx,y,Re,uei,duedsi,ueip1,duedsip1);

tsep = 0;

% 1st stage

yt(1) = y(1);
yt(2) = y(2);

H1 = yt(2);
H = HofH1(H1);

if H == 3 
   tsep = 1;
end;

if tsep ==0

  rtheta = Re*uei*yt(1);

  yp(1) = -(H+2)*yt(1)*duedsi./uei + 0.5*cfturb(rtheta,H);
  yp(2) = -H1*(duedsi./uei + yp(1)./yt(1))+0.0306*(H1-3).^(-0.6169)./yt(1);

  yt(1) = y(1) + dx*yp(1);
  yt(2) = y(2) + dx*yp(2);

  ys(1) = y(1) + 0.5*dx*yp(1);
  ys(2) = y(2) + 0.5*dx*yp(2);

  % 2nd stage

  H1 = yt(2);
  H = HofH1(H1);

  if H == 3 
    tsep = 1;
  end;

  if tsep ==0

     rtheta = Re*ueip1*yt(1);

     yp(1,1) = -(H+2)*yt(1)*duedsip1./ueip1 + 0.5*cfturb(rtheta,H);
     yp(2,1) = -H1*(duedsip1./ueip1 + yp(1)./yt(1))+0.0306*(H1-3).^(-0.6169)./yt(1);

     ynp1(1) = ys(1) + 0.5*dx*yp(1);
     ynp1(2) = ys(2) + 0.5*dx*yp(2);

  end;

end;

if tsep ==1
  ynp1(1) = -2; % so that H == 3 in the main
  ynp1(2) = 0;
end;

return

%%%% FUNCTION %%%%

function cf = cfturb(rtheta,H);

cf = 0.246*(10.^(-0.678*H))*rtheta.^(-0.268);

return

%%%% FUNCTION %%%%

function cd = sy(thup,Hup,ueteup,thlo,Hlo,uetelo);

cd = 2*thup*(ueteup).^((Hup+5)/2) + 2*thlo*(uetelo).^((Hlo+5)/2);

cd = floor(round(10000*cd))/10000;

return

%%%% FUNCTION %%%%

function z  = naca4(x,Extra);

%CREATES NACA4 AIRFOILS FROM INPUT

% Decoding arguments
eps = x(1)./100;
p = x(2)./10;
to = x(3)./100;

nbpo2 = Extra(1);
c = 1;

% x distribution   
beta = (0:(pi./nbpo2):pi).'; 
xc = c.*(1-.5.*(1-cos(beta)));   

% Thickness distribution = f(x)
thdis=5.*to.*c.*(0.2969.*sqrt(xc./c)-0.126.*xc./c-0.3537.*(xc./c).^2 +0.2843.*(xc./c).^3-0.1015.*(xc./c).^4);

% Camberline = f(x)
if p~=0 & eps~=0 
 I1=find(xc(1:nbpo2+1)<=p*c);
 I2=find(xc(1:nbpo2+1)>p*c);
 camberline(I1,1) = (eps.*xc(I1))./(p.^2).*(2.*p-xc(I1)./c);
 camberline(I2,1) = (eps.*(c-xc(I2)))./(1-p).^2.*(1+xc(I2)./c-2.*p);
end;

% Airfoil = camberline and thickness
if p==0 | eps==0
    xupper =  xc;
    yupper =  thdis;
    xlower =  xc;
    ylower = -thdis;
else
    theta(I1,1) = atan(2.*eps./p.*(-xc(I1)./(c.*p)+1));
    theta(I2,1) = atan(2.*eps./(1-p.^2).*(p-(xc(I2)./c)));
    xupper = xc         - thdis.*sin(theta);
    yupper = camberline + thdis.*cos(theta);
    xlower = xc         + thdis.*sin(theta);
    ylower = camberline - thdis.*cos(theta); 
end;

z = [xupper yupper ; xlower(nbpo2:-1:1,1) ylower(nbpo2:-1:1,1)];

return