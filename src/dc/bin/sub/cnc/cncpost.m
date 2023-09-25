function cncpost(fid,format,X,Y,Z,data,varargin),
%  CNCPOST - Output CNC data to file
%    CNCPOST(FID,FORMAT,X,Y,Z,DATA) creates formatted GCODE for the
%    toolpath given in vectors X, Y, and Z. FID is an integer file
%    identifier obtained from fopen.  Optionally, FID can be set to 1 for
%    printing to the screen. DATA is a structure containing 4 fields:
%    filename, rpm, feedrate, and units.  
%
%    Example:
%      load samplePath
%      fid = fopen(filename,'w');
%      cncpost(fid,'X%1.6g Y%1.6g Z%1.6g\n',X,Y,Z,Data)
%      cncpost(1,'X%1.6g Y%1.6g Z%1.6g\n',X,Y,Z,Data)
%      fclose(fid)
%      open('test.nc')

% GCODE REFERENCE
% G00 (Rapid Traverse)
% G01 (Linear Interpolation)
% G40 (Cancel cutter compensation)
% G80 (Cancel Drilling Cycle)
% G90 (make parameters absolute)
% G70 (make parameters inch units)
% G71 (make parameters mm units)
% G53 (set to machine coordinates)

if nargin < 6,
    data.filename = 'test.nc';
    data.rpm = 1800;
    data.feedrate = 30;
    data.units = 'in';
end

filename = data.filename;
rpm      = data.rpm;
feedrate = data.feedrate;
units    = data.units;

switch lower(units),
    case 'in',
        units = 0;
        zoff = .5;
    case 'mm',
        units = 1;
        zoff = 10;
    case 'm',
        units = 1;
        X = X*1000;
        Y = Y*1000;
        Z = Z*1000;
        zoff = 10;
    otherwise,
        units = 0;
        zoff = .5;
end
zoff = (max(Z) + zoff);

start = 1:2;
X = X(:); Xstart = X(start); X(start) = [];
Y = Y(:); Ystart = Y(start); Y(start) = [];
Z = Z(:); Zstart = Z(start); Z(start) = [];

G = ones(size(X));

% Set repeated vector components to NaN
X([1; diff(X)] == 0) = nan;
Y([1; diff(Y)] == 0) = nan;
Z([1; diff(Z)] == 0) = nan;

% Remove repeated vectors
rm = sum(isnan([X Y Z]),2) == 3;
X(rm) = [];
Y(rm) = [];
Z(rm) = [];
G(rm) = [];

% Line numbering (Start renumbering after 9999)
startat = 7;
N = [startat:length(X)+startat]';
N = mod(N,10000);
Nend = N(end);
N = N(1:end-1);

% Default format
if isempty(format),
    format = 'N%04g G%02.0f X%0.4f Y%0.4f Z%0.4f\n';
end

% Prefix with line numbering and GCODE
format = ['N%04g G%02.0f ' format];
values = [N G X Y Z]';

% Toolpath Format: N0000 G01 X0.000000 Y0.000000 Z0.000000
str = sprintf(format,values);

% Remove repeat vector components to save on file size
str = strrep(str,'XNaN ','');
str = strrep(str,'YNaN ','');
str = strrep(str,'ZNaN','');
path = str;

% Header
date = datestr(now,'dd-mm-yy');
time = datestr(now,'HH:MM');
header  = sprintf('( %s DATE=%s  TIME=%s )\n',upper(filename),date,time);
header2 = sprintf('( Toolpath created using cncpost by: Daniel Claxton )\n\n');
header = [header header2];

Pstart = [Xstart(2) Ystart(2) Zstart(2)];

pre1   = sprintf('N0000 G7%1.0f\n',units);                         % Seet Units G70(in) or G71(mm)
pre2   = sprintf('N0001 G00 G40 G80 G90\n');                       % Intial G codes
pre3   = sprintf('N0002 T1\n');                                    % Select Tool 1
rpm    = sprintf('N0003 G00 Z%0.4f S%1.0f M03\n',zoff,rpm);        % Raise Tool, Set Spindle Speed, Start
begin  = sprintf('N0004 G00 X%0.4f Y%0.4f\n',Xstart(1),Ystart(1)); % Go To X,Y
begin2 = sprintf('N0005 G01 Z%0.4f E%1.0f\n',Zstart(1),feedrate(2));        % Go To Z, set Z feedrate
feed   = sprintf('N0006 G01 X%0.4f Y%0.4f Z%0.4f F%1.0f\n',Pstart,feedrate(1)); % Set XY feedrate
finish = sprintf('N%04g M30\n',Nend);                              % Turn off

% Concatenate all strings
pre = [pre1 pre2 pre3];
begin = [begin begin2];
out = [header pre rpm begin feed path finish];

fprintf(fid,'%s',out);

