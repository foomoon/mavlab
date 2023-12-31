function nrb2iges(obj,filename,varargin)

[obj,Entity,Num] = getEntities(obj);


% Decompose file name
[path file ext]=fileparts(filename);

% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
%                           Start (S-SECTION)
% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
S = ['Matlab to IGES converter. Written by Daniel Claxton',...
      blanks(21) 'S0000001'];
% SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS





% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
%                        Setup Header (G-SECTION)
% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
G.paramdelim     = ',';                           % Parameter Deliminator Character 
G.recorddelim    = ';';                           % Record Delimiter Character
G.sendid         = 'Sender ID';                   % Product ID from Sender
G.receiveid      = 'Receiver ID';                 % Product ID for Receiver
G.filename       = [file '.igs'];                 % File Name
G.systemid       = computer;                      % System ID
G.preprocessor   = 'Matlab -> IGES';              % Pre-processor Version
G.intbits        = 16;                            % Number of Bits for Integers
G.precision1     = 06;                            % Single Precision Magnitude
G.significant1   = 15;                            % Single Precision Significance
G.precision2     = 13;                            % Double Precision Magnitude
G.significant2   = 15;                            % Double Precision Significance
G.scale          = 1.0;                           % Model Space Scale
G.unitflag       = 3;                             % Unit Flag (1 = inches, 3 = look to units)
G.units          = 'IN';                          % Model Units 'MM' or 'IN' (Default)
G.lineweights    = 8;                             % Maximum number of line weights
G.linewidth      = 0.016;                         % Maximum line width
G.timestamp      = datestr(now);                  % Time stamp of creation
G.resolution     = 1E-4;                          % Minimum User-inted Resolution
G.max            = getMax(obj);                   % Approximate Maximum Coordinate
G.author         = 'D. Claxton dclaxton@ufl.edu'; % Author
G.company        = 'University of Florida';       % Author's Organization
G.version        = 11;                            % - USPRO/IPO-100 (IGES 5.2) [USPRO93]';  % IGES Version Number  ** prob not right **
G.standard       = 3;                             % - ANSI'; % Drafting Standard Code
G.date           = datestr(today);                % Model Creation/Change Date

% Check inputs
if nargin < 3
    GG = G;
else
    % Set above parameters from input parameter value pairs
    GG = setParams(G,varargin{:});
end

G = G2str(GG);
% GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG






%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
%                              P-SECTION
%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
P = zeros(10,Num);
L = zeros(Num,1);
for i = 1:Num;
    switch Entity(i)
        case 128
            P(1,i) = Entity(i);             % Entity Type.  128 = Rational spline surface
            P(2,i) = obj{i}.number(1,1)-1;  % Number of points in u direction (doesn't work without -1)
            P(3,i) = obj{i}.number(1,2)-1;  % Number of points in v direction
            P(4,i) = obj{i}.order(1,1)-1;   % Degree_u  (degree = order -1)
            P(5,i) = obj{i}.order(1,2)-1;   % Degree_v
            P(6,i) = 0;                     % PROP1 Closed_u (0 = Not Closed)
            P(7,i) = 0;                     % PROP2 Closed_v (0 = Not Closed) ****change to 1 for FEA Limb****
            P(8,i) = 0;                     % PROP3 (1 = Polynomial i.e. all weights equal, 0 = rational) % *Note, we choose 0 because it is general case
            P(9,i) = 0;                     % PROP4 1st Direction periodicity (0 = Non-periodic)
            P(10,i)= 0;                     % PROP5 2nd Direction periodicity (0 = Non-periodic)
            Uknots = obj{i}.knots{1}(:);
            Vknots = obj{i}.knots{2}(:);
            Knots = [Uknots; Vknots];
        case 126            
            P(1,i) = Entity(i);             % Entity Type.  126 = Rational spline curve
            P(2,i) = obj{i}.number - 1;     % Number of points - 1
            P(3,i) = obj{i}.order  - 1;     % Degree of curve
            P(4,i) = 0;                     % PROP1 Planar curve (0 = Not Planar***
            P(5,i) = 0;                     % PROP2 Closed curve (0 = Not Closed)
            P(6,i) = 0;                     % PROP3 (1 = all weights equal, 0 = rational)***
            P(7,i) = 0;                     % PROP4 Periodic (0 = Not periodic)
            
            % ***Note, we choose 0 because it is general case, it just requires us to
            % write more information to the IGES file. But we can live with that to
            % simplify things.
            Knots = obj{i}.knots(:);
        otherwise
    end

    if Entity(i) == 126 || Entity(i) == 128
        X = squeeze(obj{i}.coefs(1,:,:));   % Extract X Data
        Y = squeeze(obj{i}.coefs(2,:,:));   % Extract Y Data
        Z = squeeze(obj{i}.coefs(3,:,:));   % Extract Z Data
        W = squeeze(obj{i}.coefs(4,:,:));   % Extract Weights from coefficients

        Coefs = [X(:) Y(:) Z(:)]';
        Weights = W(:);
        Coefs = Coefs(:);

        Data{i} = [Knots; Weights; Coefs];
    end


end
[P,L] = P2str(P,Data,GG,Num);
%PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP





% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
%                               D-SECTION
% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

D = zeros(15,Num);
for i=1:Num
    D(1,i) = Entity(i);                 % Entity Type.
    D(2,i) = 1+L(i);                    % Data Start line
    D(3,i) = 0;                         % Structure
    D(4,i) = 1;                         % Line Font Pattern (1 = Solid)
    D(5,i) = 0;                         % Level
    D(6,i) = 0;                         % View
    D(7,i) = 0;                         % Transformation Matrix
    D(8,i) = 0;                         % Label Display
    D(9,i) = 0;                         % Blank Status (0 = Visible)
    D(10,i)= 0;                         % Subord. Entity Switch (0 = Independant)
    D(11,i)= 0;                         % Entity Use Flag (0 = Geometry)
    D(12,i)= 1;                         % Hierarchy ( 1 = Global defer)
    D(13,i)= 0;                         % Line Weight Number
    D(14,i)= L(i+1);                    % Data end line (Will be set later)
    D(15,i)= 0;                         % Form Number (9 = General Quadratic Surface), 0 = none of above (1-9) options
end
D = D2str(D);
% DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD






% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
%                           T - SECTION
% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
T = [size(S,1),size(G,1),size(D,1),size(P,1)];
T = sprintf('S%07.0fG%07.0fD%07.0fP%07.0f%sT%07.0f',T,blanks(40),1);
% TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

chk = [size(S); size(G); size(D); size(P); size(T)];

if sum(diff(chk(:,2))) ~= 0
   error(['Something went wrong writing the IGES file: ' file '.igs']) 
end

% Write to file
str = [S; G; D; P; T];

fid = fopen([path file '.igs'],'w');


for i=1:size(str,1)-1
    fprintf(fid,'%s\n',str(i,:));
end
fprintf(fid,'%s',str(i+1,:));

fclose(fid);













function [obj,ent,n] = getEntities(obj)
if ~iscell(obj)
    obj = {obj};
end
n = length(obj);
for i=1:n
    if isstruct(obj{i})
        if iscell(obj{i}.knots)
            ent(i) = 128;
        else
            ent(i) = 126;
        end
    else
        [rows,cols] = size(obj{i});
        if rows == 1
            ent(i) = 116;
        elseif cols == 1
            ent(i) = 116;
            obj{i} = obj{i}';
        else
            ent(i) = 110;
            if cols > rows
                obj{i} = obj{i}';
            end
        end
    end
end




function G = setParams(G,p)
n = length(p);
n = n/2;
p = reshape(p,n,2);
for i=1:n
    G = setfield(lower(G,p{i,1}),p{i,2});
end




function s = wrap(str,delim,n)
i = 1;
while true
    if length(str) > n
        tmp = str(1:n);
        d = findstr(tmp,delim);
        d = d(end);
        if d == length(tmp)
            str = str(n+1:end);
        else
            str = [tmp(d+1:end) str(n+1:end)];
        end
        tmp = tmp(1:d);
        tmp(d+1:n) = blanks(n-d);
        s(i,:) = tmp;
        i = i+1;
    else
        m = length(str);
        stop = n - m;
        tmp = [str blanks(stop)];
        s(i,:) = tmp;
        break
    end
end


function str = D2str(D)

k = 1;
[m,n]=size(D);
for Nent=1:n
    Entity = D(1,Nent);
    %     switch Entity
    %         case 128
    %         case 126
    %         case 110
    %         case 116
    %         otherwise
    %     end
    tmp1 = sprintf('%8.0f', D(1:7,Nent));
    tmp2 = sprintf('%8.0f',[D(1,Nent); D(12:15,Nent)]);
    str(k,:)   = [tmp1 blanks(23-8-7)];
    k = k + 1;
    str(k,:) = [tmp2 blanks(39-8-7)];
    k = k + 1;
end

[m,n] = size(str);
on = true;
for i=1:m
    if on
        count(i,:) = sprintf('% 8.0fD%07.0f',i,i);
    else
        count(i,:) = sprintf('% 8.0fD%07.0f',0,i);
    end
    on = ~on;
end

str = [str count];


function [str,L] = P2str(P,Data,G,Nentities)
pd = G.paramdelim;
rd = G.recorddelim;
str = [];
L(1) = 0;
for Nent=1:Nentities
    switch P(1,Nent)
        case 128
            plen = 1:10;
        case 126
            plen = 1:7;
        case 110
        case 116
        otherwise
    end
    str1 = sprintf(['%1.0f' pd],P(plen,Nent));
    str2 = sprintf(['%1.6f' pd],Data{Nent});
    str3 = ['0., 1., 0., 1.' rd];

    tmp = [str1 str2 str3];
    nsp = 7;

    % Wrap String to be
    tmp = wrap(tmp,pd,71-nsp);

    [m,n] = size(tmp);
    L(Nent+1,1) = m;

    for i=1:n-64
        tmp(:,end+1) = ' ';
    end
    
    for j=1:m
        count(j,:) = sprintf('% 8.0fP%07.0f',2*Nent-1,j+L(Nent));
    end

    tmp = [tmp count];
    
    str = [str; tmp];
    
    count = []; % Clear count
   
end







function str = G2str(G)
fnames = fieldnames(G);
n = length(fnames);
str = cell(20,1);
pd = G.paramdelim;
rd = G.recorddelim;

k = 1;
for i=1:n
   tmp = getfield(G,fnames{i});
   if isstr(tmp)
       tmp = sprintf('%1.0fH%s%s',length(tmp),tmp,pd);
   else
       if tmp-fix(tmp)
           tmp = sprintf('%1.6f%s',tmp,pd);
       else
           tmp = sprintf('%1.0f%s',tmp,pd);
       end
   end
   if length(str{k}) + length(tmp) > 79-15
       k = k + 1;
   end
   
   str{k,1} = [str{k} tmp];
   if i == n
       str{k}(end) = rd;
   end
end

last = 0;
while true
    if isempty(str{last+1})
        break
    end
    last = last+1;
end

str = str(1:last,1);
n = length(str{end,1});
stop = 80-8;
str{end}(n+1:stop) = blanks(stop-n);
str = char(str);

for i=1:last
    count(i,:) = sprintf('G%07.0f',i);
end

str = [str count];



function m = getMax(obj)

for i=1:length(obj)
    if isstruct(obj{i})
        m(i) = max(obj{i}.coefs(:));
    else
        m(i) = max(obj{i}(:));
    end
end
m =  max(m);
