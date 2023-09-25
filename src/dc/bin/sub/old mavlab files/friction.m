function cf = friction(Xme,Rex,varargin)
% FRICTION Computation Skin Friction Coeficient
%   FRICTION returns the skin friction coefficient of a surface at a given
%   mach and Reynolds number.  It handles laminar flow, turbulent flow and
%   composite flow (when laminar and turbulent both exist).
%
%   FRICTION(Mach,Re,flowtype) returns the skin friction coefficient for
%   the given flow type.  Options for "flowtype" are as follows:
%      
%       flowtype
%       --------
%       'l' - Laminar Flow
%       't' - Turbulent Flow
%        X  - Composite flow (0 < X < 1)
%
%   For the composite flow case, X represents the normalized distance from
%   the leading edge of the flow where separation occurs.  For instance, if
%   separation occurs and the middle of the chord of the wing, then X = .5.
%
%   Input:
%       Mach:       (scaler) Mach Number
%       Re:         (scaler) Reynolds Number
%       Flowtype:   (scaler/string) Flow Type
%
%   Output:
%       Cf:     (scaler) Skin Friction Coefficient
%
%   *Note: This function is valid for Mach <= 3.
%
%   See also FORMFACTOR

% This code is modified from an existing fortran code written by:
% W. Mason,   October 21, 1989
% Department of Aerospace and Ocean Engineering
% Virginia Tech, Blacksburg, VA 24061
% mason@aoe.vt.educ     
%
% Adapted by: 
% Daniel Claxton,   February 26, 2005
% Department of Mechanical & Aerospace Engineering
% University of Florida, Gainesville, FL 32611
% dclaxton@ufl.edu

% See http://www.aoe.vt.edu/~mason/Mason_f/MRsoft.html#SkinFriction
% For Fortran Code and Manual for it

TwTaw = 1;   % Adiabatic Wall Temperature Ratio (Assume = 1 for M<3)

cfl_x = lamcf(Xme,Rex,TwTaw);
cft_x = turbcf(Xme,Rex,TwTaw);
if Xme > 3,
    disp(' ')
    disp(' Warning: Results may not be acurate for Mach numbers greater than 3');
end

if strcmp(lower(varargin{1}),'l'),
    cf = cfl_x;
elseif strcmp(lower(varargin{1}),'t'),
    cf = cft_x;
elseif ~isstr(varargin{1}),
    sepx = varargin{1};
    Rec = sepx*Rex;
    cfl_c = lamcf(Xme,Rec,TwTaw);
    cft_c = turbcf(Xme,Rec,TwTaw);
    cfc = cft_x - (sepx)*(cft_c - cfl_c);
    cf = cfc;
else
    cf = -1;
end


function CF_lam = lamcf(Xme,Rex,TwTaw)
G     = 1.4;        % Specific Heat Ratio
Pr    = .72;        % Prandl Number
R     = sqrt(Pr);   % Recovery Factor
Te    = 390;        % Edge Temperature
Tk    = 200;        % 

TwTe  = TwTaw*(1 + R*((G - 1)/2)*Xme^2);
TstTe = 0.5 + 0.039*Xme^2 + 0.5*TwTe;
Cstar = sqrt(TstTe)*((1 + (Tk/Te))/(TstTe + (Tk/Te)));
CF_lam  = 2*.664*sqrt(Cstar)/sqrt(Rex);



function CF_turb = turbcf(Xme,Rex,TwTaw)
epsmax = 0.1e-8;
G      = 1.4;
r      = 0.88;
Te     = 222.0;

xm    = (G - 1)/2*Xme^2;
TawTe = 1. + r*xm;
F     = TwTaw*TawTe;
Tw    = F * Te;
A     = sqrt(r*xm/F);
B     = (1. + r*xm - F)/F;
denom = sqrt(4*A^2 + B^2);
Alpha = (2.*A^2 - B)/denom;
Beta  = B/denom;
Fc    = ((1.0 + sqrt(F))/2.0)^2;

if Xme >= .1,
    Fc = r*xm/(asin(Alpha) + asin(Beta))^2;
end

Xnum   = (1. + 122./Tw*10^(-5/Tw));
Denom  = (1. + 122./Te*10^(-5/Te));
Ftheta = sqrt(1./F)*(Xnum/Denom);
Fx     = Ftheta/Fc;

RexBar = Fx * Rex;

Cfb    = 0.074/(RexBar^0.20);

iter = 1;
epsi =1;
while (epsi > epsmax),
    iter = iter + 1;
    if iter > 200, break; end
    Cfo    = Cfb;
    Xnum   = 0.242 - sqrt(Cfb)*log10(RexBar*Cfb);
    Denom  = 0.121 + sqrt(Cfb)/log(10);
    Cfb    = Cfb*(1.0 + Xnum/Denom);
    epsi    = abs(Cfb - Cfo);
end

CF_turb     = Cfb/Fc;


