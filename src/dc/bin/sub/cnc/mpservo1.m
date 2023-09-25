function mpservo1(X,Y,Z,fullpath,tool,diameter,feedrate)
% MPSERVO1  Toolpath post-processor.
%   MPSERVO1(X,Y,Z,filename,tool,diameter,feedrate) outputs  a
%   toolpath (.NC) file for the toolpath given in vectors X,Y,Z. Filename
%   is the path and name of the file it will be saved to.  If it does not
%   exist, it will be created.  Tool is a string specifying the type of
%   tool used to cut, ie '1/2 BALL ENDMILL'.  Diameter specifies the
%   diameter of the cutting tool.  Feedrate specifies the rate at which the
%   machine will cut in inches/minute. 
%
%   MPSERVO1 returns X,Y,Z coordinates as a set of column vectors.
%   
%   See also CNC, MAVLAB, FLOWLINE


% Round Numbers to Precision of Machine
n = 4;                                  % Precision
X = round(X*10^n)/(10^n);
Y = round(Y*10^n)/(10^n);
Z = round(Z*10^n)/(10^n);

% Header Data
date = datestr(today,20);
date(findstr(date,'/')) = '-';
time = datestr(now,15);

[p,filename]=fileparts(fullpath);

% ---- Start Writing File ----
fid = fopen([fullpath '.nc'],'w');

% Header
fprintf(fid,'( %s NC   DATE=DD-MM-YY - %s  TIME=HH:MM - %s )\n',upper(filename),date,time);
fprintf(fid,'(  %s   TOOL - 1 DIA. OFF. - 21 LEN. - 2 DIA. - %0.1f )\n',tool,diameter);

% Start
fprintf(fid,'N100 G70\n');
fprintf(fid,'N102 G00 G40 G80 G90\n');
fprintf(fid,'N104 T1\n');
fprintf(fid,'N106 G00 G90 X%0.4f Y%0.4f S1069 M03\n',X(1),Y(1));
fprintf(fid,'N112 G01 Z%0.4f E%1.2f\n',Z(1),feedrate);
fprintf(fid,'N114 X%0.4f Y%0.4f Z%0.4f F%1.2f\n',X(2),Y(2),Z(2),feedrate);

count = 116;
for i=3:length(X),
    if X(i-1) ~= X(i) & Y(i-1) ~= Y(i) & Z(i-1) ~= Z(i),
        % If X,Y,Z are unique with respect to the last set, print each of them
        fprintf(fid,'N%1.0f X%0.4f Y%0.4f Z%0.4f\n',count,X(i),Y(i),Z(i));
    elseif X(i-1) == X(i) & Y(i-1) == Y(i) & Z(i-1) ~= Z(i),
        % If X & Y are same as last coordinate don't repeat
        fprintf(fid,'N%1.0f Z%0.4f\n',count,Z(i));
    elseif X(i-1) == X(i) & Y(i-1) ~= Y(i) & Z(i-1) == Z(i),
        % If X & Z are same as last coordinate don't repeat
        fprintf(fid,'N%1.0f Y%0.4f\n',count,Y(i));
    elseif X(i-1) ~= X(i) & Z(i-1) == Z(i) & Y(i-1) == Y(i),
        % If Y & Z are same as last coordinate don't repeat
        fprintf(fid,'N%1.0f X%0.4f\n',count,X(i));
    elseif X(i-1) == X(i) & Y(i-1) == Y(i) & Z(i-1) == Z(i),
        % If X,Y,Z are all the same as last coordinate
        % Do nothing
        count = count - 2;          % null out the count for this iteration
    elseif X(i-1) == X(i),
        % If X is the same as last coordinate don't repeat
        fprintf(fid,'N%1.0f Y%0.4f Z%0.4f\n',count,Y(i),Z(i));
    elseif Y(i-1) == Y(i),
        % If Y is the same as last coordinate don't repeat
        fprintf(fid,'N%1.0f X%0.4f Z%0.4f\n',count,X(i),Z(i));
    elseif Z(i-1) == Z(i),
        % If Z is the same as last coordinate don't repeat
        fprintf(fid,'N%1.0f X%0.4f Y%0.4f\n',count,X(i),Y(i));    
    end
    count = count+2;
    if count == 10000, count = 100; end
end

fprintf(fid,'N%1.0f G00 Z.0472\n',count); count = count+2;
if count == 10000, count = 100; end
fprintf(fid,'N%1.0f Z2.5\n',count); count = count+2;
if count == 10000, count = 100; end
fprintf(fid,'N%1.0f M05\n',count); count = count+2;
if count == 10000, count = 100; end
fprintf(fid,'N%1.0f G91 Z0.\n',count); count = count+2;
if count == 10000, count = 100; end
fprintf(fid,'N%1.0f M30\n',count);

fclose(fid);
% ---- Finish Writing ----

