function mat2avl(fid,X,Y,Z,varargin),


% Shift Wing in the chord-wise direction such that the centerleading edge 
% is at [0,0,0]
X = X - min(X(:));

name = varargin{1};
mach = varargin{2};
Npan = varargin{3};
Nspan = Npan(1);
Nchord = Npan(2);

% % REFORMAT DATA TO BE CONSISTENT WITH AVL
% XX = Y(21:end,1);
% YY = X(21:end,1);
% ZZ = Z(21:end,1);
% 
% % CALCULATE ROOT CHORDS
% chord = XX - Y(21:end,end);

% EXTRACT DATA FROM WING
[XX YY ZZ CC AA] = extract(X,Y,Z);
% XX = -(XX-max(XX));

% CALCULATE REFERENCE SURFACE AREA, POSITION, SPAN & CHORD
[Sref,Xref,Bref,Cref] = getRef(X,Y,Z);
% Xref = 0;

% FORMAT AIRFOILS
% [foilX,foilY] = blendfoils(wing);
[foilX,foilY] = getFoils(X,Y,Z,AA,CC);
foilX = [flipud(foilX); foilX(2:end,:)];
foilY = [flipud(foilY); foilY(2:end,:)];

% -----------------
% START HEADER DATA
% -----------------
% The input file begins with the following information in the first 5 non-blank,
% non-comment lines:
% 
% Abc...              | case title
% 
% #                   | comment line begins with "#" or "!"
% 
% 0.0                 | Mach
% 1     0     0.0     | iYsym  iZsym  Zsym
% 4.0   0.4   0.1     | Sref   Cref   Bref
% 0.1   0.0   0.0     | Xref   Yref   Zref
% 0.020               | CDp  (optional)

SYM = [0 0 0];
SCB = [Sref Cref Bref];
XYZ = [Xref 0 0];
% ref = [SYM; SCB; XYZ];

fprintf(fid,'%s\n\n',name);
fprintf(fid,'#Mach\n\n');
fprintf(fid,'%1.4f\n',mach);
fprintf(fid,'#IYsym   IZsym   Zsym\n');
fprintf(fid,'%1.4f %1.4f %1.4f\n',SYM);
fprintf(fid,'#Sref    Cref    Bref\n');
fprintf(fid,'%1.4f %1.4f %1.4f\n',SCB);
fprintf(fid,'#Xref    Yref    Zref\n');
fprintf(fid,'%1.4f %1.4f %1.4f\n',XYZ);
% fprintf(fid,'%s',mat2str2(ref,'%1.3f\t'));
if size(varargin,2) > 3,
    CDp  = varargin{4}; %Profile Drag
    fprintf(fid,'#Profile Drag\n');
    fprintf(fid,'%1.5f\n',CDp);
end

% -----------------
% END HEADER DATA
% -----------------



% -----------------
% PARAMETERS
% -----------------
Nsect  = length(XX);
% Nchord = 8;
Cspace = 1;
% Nspan  = 1;
Sspace = 1;

% SCALE
Xscale = 1;
Yscale = 1;
Zscale = 1;
Scale  = [Xscale Yscale Zscale];

% TRANSLATE
DX = 0;
DY = 0;
DZ = 0;
Translate = [DX DY DZ];

% GLOBAL ANGLE OF INCIDENCE
Alpha = 0;

% -----------------
% PARAMETERS
% -----------------

fprintf(fid,'\n#====================================================================\n');

% -------------
% START SURFACE
% -------------

% DEFINE SURFACE PARAMETERS
fprintf(fid,'SURFACE\n');
fprintf(fid,'Main Wing\n');
fprintf(fid,'#Nchordwise  Cspace   Nspanwise   Sspace\n');
fprintf(fid,'%1.3f        %1.3f    %1.3f       %1.3f\n',[Nchord Cspace Nspan Sspace]);
fprintf(fid,'#\n');
% % DEFINE INDEX PARAMETERS
% fprintf(fid,'INDEX\n');
% fprintf(fid,'1\n');

% DEFINE YDUPLICATE
fprintf(fid,'YDUPLICATE\n');
fprintf(fid,'%1.0f\n',0);
fprintf(fid,'#\n');

% DEFINE SCALE
fprintf(fid,'SCALE\n');
fprintf(fid,'#X     Y      Z\n');
fprintf(fid,'%1.4f %1.4f %1.4f\n',Scale);
fprintf(fid,'#\n');

% DEFINE TRANSLATE
fprintf(fid,'TRANSLATE\n');
fprintf(fid,'#X     Y      Z\n');
fprintf(fid,'%1.4f %1.4f %1.4f\n',Translate);
fprintf(fid,'#\n');

% DEFINE ANGLE
fprintf(fid,'ANGLE\n');
fprintf(fid,'%1.2f\n',Alpha);
fprintf(fid,'\n');

% -------------
% END SURFACE
% -------------

% -------------------------
% LOOP THROUGH EACH SECTION
% -------------------------
if ~isdir([cd filesep name filesep 'AIRFOILS'])
    mkdir([name '\AIRFOILS']);
end
for i=1:Nsect,
    % START SECTION
    fprintf(fid,'#-------------------------------------------------------------\n');
    fprintf(fid,'SECTION\n');
%     fprintf(fid,'%s',mat2str2([XX(i) YY(i) ZZ(i) CC(i) AA(i) Nspan Sspace],'%1.4f '));
    fprintf(fid,'#Xle    Yle     Zle     Chord   Ainc    Nspanwise Sspace\n');
%     fprintf(fid,'%1.4f  %1.4f  %1.4f  %1.4f  %1.4f  %1.4f    %1.4f\n',[XX(i) YY(i) ZZ(i) CC(i) AA(i) 0*Nspan 0*Sspace]);
    fprintf(fid,'%1.4f  %1.4f  %1.4f  %1.4f  %1.4f\n',[XX(i) YY(i) ZZ(i) CC(i) AA(i)]);
    
    % AIRFOIL ( Doesn't seem to work, so we have to use AFILE )
%      fprintf(fid,'AIRFOIL\n');
%      fprintf(fid,'%1.4f %1.4f\n',[foilX(:,i),foilY(:,i)]');

    % AFILE
    file = [name '\AIRFOILS\' num2str(i) '.dat'];
    printairfoil(foilX(:,i),foilY(:,i),file,name);
    fprintf(fid,'AFILE\n');
    fprintf(fid,'%s\n',file);
    fprintf(fid,'\n');
end



function printairfoil(X,Y,file,name),
fid = fopen(file,'w');
fprintf(fid,'%s\n',file);
% fprintf(fid,'%s',mat2str2([X,Y],'%1.5f\t'));
fprintf(fid,'%1.5f %1.5f\n',[X(:) Y(:)]');
fclose(fid);


function [xl yl zl c alpha] = extract(X,Y,Z)
xl = X(:,1);
yl = Y(:,1);
zl = Z(:,1);
xt = X(:,end);
yt = Y(:,end);
zt = Z(:,end);
c = sqrt((xt-xl).^2 + (yt-yl).^2 + (zt-zl).^2);
h = zt-zl;
alpha = -(180/pi)*tan(h./c);
% for i= 2:length(alpha)
%     alpha(i) = alpha(i) + alpha(i-1);
% end


function [foilX,foilY] = getFoils(X,Y,Z,AA,CC)
for i=1:size(X,1)
    theta = (pi/180)*AA(i);
    T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    x = (X(i,:)-X(i,1));  
    y = (Z(i,:)-Z(i,1));  
    n = T*[x; y];
    foilX(:,i) = n(1,:)'/CC(i);
    foilY(:,i) = n(2,:)'/CC(i);
end



function [Sref,Xref,Bref,Cref] = getRef(X,Y,Z)
% Sref = 2*polyarea([Y(:,1); Y(end:-1:1,end)],[X(:,1); X(end:-1:1,end)]);
% [meanAC Cref] = mac(Y(:,1),X(:,1),Y(:,end),X(:,end));
% Xref = (X(round(end/2),1) - meanAC);
% Bref = 2*(max(Y(:))-min(Y(:)));
% Cref = abs(Cref);
X = -X;
XX = Y;
Y = X;
X = XX;
Y = Y-max(Y(:));
% Get Border Data 
x1 = fliplr(X(1,:)); x2 =  X(:,1)'; x3 = X(end,:); x4 = flipud(X(:,end))';
y1 = fliplr(Y(1,:)); y2 =  Y(:,1)'; y3 = Y(end,:); y4 = flipud(Y(:,end))';
x = [x1 x2 x3 x4];
y = [y1 y2 y3 y4];
xl = X(:,1);
yl = Y(:,1);
xt = X(:,end);
yt = Y(:,end);

Area = 2*polyarea(x,y);
[mAC mChord] = mac(xl,yl,xt,yt,Area);
% AC = abs(Y(1,1) - mAC)/wing.chord;
% AR = wing.span/mChord; %(wing.span)^2/Area;

Sref = Area;
Cref = mChord;
Xref = Y(1,1);
Bref = 2*(max(X(:))-min(X(:)));


    
        
function [ac,mean_ac] = mac(xl,yl,xt,yt,A)

% figure; plot(xl,yl,'r',xt,yt,'g')
x = linspace(xl(1),xl(end),5000)';
yl = interp1(xl,yl,x);
yt = interp1(xt,yt,x);
xl = x;

c = yl-yt;
dx = [0; diff(xl)];
% mean_ac = sum(c.^2)/sum(c);
mean_ac = 2*sum(dx.*c.^2)/A;
xb = sum(yl.*c)/sum(c);
ac = xb - (mean_ac/4);





% function [X Y Z C A] = extract2(wing),
% 
% 
% % Make Data Names Easier to Use (Because I'm Lazy)
% span = wing.span;
% chord = wing.chord;
% theta = wing.dihedral;
% polyhedral = wing.polyhedral;
% sweep_angle = wing.sweep;
% twist_angle = wing.twist;
% curvy = wing.curvature;
% curvetype = wing.curvetype;
% 
% % Don't ask why I stored the planform data this way
% xl=wing.ledge(:,1);
% yl=wing.ledge(:,2);
% xt=wing.tedge(:,1);
% yt=wing.tedge(:,2);
% 
% % Scale Normalized Wing Data to Correct Span & Chord
% xl = xl*span/2;
% xt = xt*span/2;
% yl = yl*chord;
% yt = yt*chord;
% 
% % Sweep Planform to Specified Angle
% [xl,yl] = sweep(xl,yl,sweep_angle);
% [xt,yt] = sweep(xt,yt,sweep_angle);
% 
% % Determine the Z-coordinates of planform
% zcurve = zoffset(xl,span/2,chord/2,curvy,curvetype);
% zdihedral = dihedral(xl,theta,polyhedral);
% ztotal = zdihedral+zcurve;
% 
% 
% % Sort Airfoil Data
% [x,foils] = blendfoils(wing);
% 
% if ~wing.edgeref
%     % Shift Airfoils to top edge reference
%     for i=1:length(wing.ledge),
%         chord_length(:,i) = (wing.ledge(:,2)-wing.tedge(:,2));
%     end
%     chord_length = wing.chord*chord_length';
%     fn = foils.*chord_length;
%     maxz = max(max(fn));
%     maxperfoil = max(fn);
%     maxoffset = (maxz-maxperfoil);
%     ztotal = ztotal + maxoffset;
% end
% 
% X = yl;
% Y = xl;
% Z = ztotal;
% A = (1-abs(X)/(span/2))*twist_angle; % Fix this
% C = yl-yt;






