function mat2igs(X,Y,Z,varargin),
% MAT2IGS IGES converter
%   MAT2IGS converts surface data into IGES format.
%
%   MAT2IGS(X,Y,Z,file) turns the the surface defined by matrices X,Y,Z
%   into line segment cross-sections defined by the rows of each matrix.
%   If file = 1, the data is printed to the screen, rather than a file.
%
%   MAT2IGS(X,Y,Z), without a filename prompts the user to select a file to
%   save to.  
%
%   Input:
%       X:      (2D - Array) x coordinates
%       Y:      (2D - Array) y coordinates
%       Z:      (2D - Array) z coordinates
%       file:   (string) save filename
%
%   Output:
%       none
%
% See also: MAVLAB, FPRINTF

if nargin < 4,
    [filename pathname] = uiputfile('*.igs','Export');
    fid=fopen([pathname filename '.igs'],'w');
else
    fid=fopen([varargin{1} '.igs'],'w');
end

% Shift everything so it's positive
X = X + abs(min(min(X)));
Y = Y + abs(min(min(Y)));
Z = Z + abs(min(min(Z)));

% Get Rid of Tips
X = X(2:end-1,:);
Y = Y(2:end-1,:);
Z = Z(2:end-1,:);

% Get Size of mesh
[m,n] = size(X);

igesline(fid,X,Y,Z);

fclose(fid);

% --- Line Entity
function igesline(fid,X,Y,Z),

[m,n] = size(X);

% S - Section
fprintf(fid,'%sS%s1\r\n',spaces(72),spaces(6));  % Nothing for now

% G - Section
fprintf(fid,'%sG%s1\r\n',spaces(72),spaces(6));  % Nothing for now

% D - Section 
entity = '110';  % Line Entity

entity_num = 1;
k = 7 - length(num2str(entity_num));
k2 = 7 - length(num2str(entity_num+1));
objn = 1;
j1 = 8 - length(num2str(objn));
fprintf(fid,'%s%s%s%1.0f%s0%s1%s0%s0%s0%s000000000D%s%1.0f\r\n',spaces(5),entity,spaces(j1),objn,spaces(7),spaces(7),spaces(7),spaces(7),spaces(7),spaces(7),spaces(k),entity_num);
fprintf(fid,'%s%s%s0%s8%s%1.0f%s0%s %s %s         D%s%1.0f\r\n',spaces(5),entity,spaces(7),spaces(7),spaces(7),2,spaces(7),spaces(7),spaces(7),spaces(7),spaces(k2),entity_num+1);

% Name
nameline = 3;
fprintf(fid,'     406%s%1.0f       0       0       0       0       0       000000000D      3\r\n',spaces(7),nameline);
fprintf(fid,'     406       0       0       1       3                                D      4\r\n');

line_num = nameline;
objn = 2;
for i=1:m,  
    for j=1:n-1,
        if i==1 & j==1,            
        else
            line_num = line_num + 2;%m*(i-1)+2*j+1;
            k = 7 - length(num2str(line_num));
            k2 = 7 - length(num2str(line_num+1));
            objn = objn + 2;%(m-1)*(i-1) + 2*j;
            j1 = 8 - length(num2str(objn));
            j2 = 8 - length(num2str(2));
            fprintf(fid,'%s%s%s%1.0f%s0%s1%s0%s0%s0%s000000000D%s%1.0f\r\n',spaces(5),entity,spaces(j1),objn,spaces(7),spaces(7),spaces(7),spaces(7),spaces(7),spaces(7),spaces(k),line_num);
            fprintf(fid,'%s%s%s0%s8%s%1.0f%s0%s %s %s         D%s%1.0f\r\n',spaces(5),entity,spaces(7),spaces(7),spaces(j2),2,spaces(7),spaces(7),spaces(7),spaces(7),spaces(k2),line_num+1);
        end    
    end
end

% P - Section
i=1; j=1;
fprintf(fid,'%s,','110');
fprintf(fid,'%1.6E,%1.6E,%1.6E,',[X(i,j); Y(i,j); Z(i,j);]);

obj_num = (m-1)*(i-1) + 2*j - 1; 
q = length(num2str(obj_num));
num_lines = m*(i-1)+2*j-1;
p = 7 - length(num2str(num_lines));

fprintf(fid,'%s%1.0fP%s%1.0f\r\n',spaces(26-q),obj_num,spaces(p),num_lines);
fprintf(fid,'%1.6E,%1.6E,%1.6E;',[X(i,j+1); Y(i,j+1); Z(i,j+1)]);

num_lines = m*(i-1)+2*j;
p = 7 - length(num2str(num_lines));

fprintf(fid,'%s%1.0fP%s%1.0f\r\n',spaces(30-q),obj_num,spaces(p),num_lines);



k3 = 7-length(num2str(num_lines+1));
fprintf(fid,'406,2,0,1H0;%s3P%s%1.0f\r\n',spaces(59),spaces(k3),num_lines+1);

for i=1:m            
    for j=1:n-1
        if i==1 & j==1,            
        else
            fprintf(fid,'%s,','110');
            fprintf(fid,'%1.6E,%1.6E,%1.6E,',[X(i,j); Y(i,j); Z(i,j);]);            
            num_lines = num_lines + 2; 
            obj_num = num_lines + 1; 
            q = length(num2str(obj_num));
            
            p = 7 - length(num2str(num_lines));
            fprintf(fid,'%s%1.0fP%s%1.0f\r\n',spaces(26-q),obj_num,spaces(p),num_lines);
            fprintf(fid,'%1.6E,%1.6E,%1.6E;',[X(i,j+1); Y(i,j+1); Z(i,j+1)]);
            p = 7 - length(num2str(num_lines + 1));
            fprintf(fid,'%s%1.0fP%s%1.0f\r\n',spaces(30-q),obj_num,spaces(p),num_lines+1);
        end        
    end        
end

% T - Section
k = 7-length(num2str(obj_num+1));
fprintf(fid,'S%s1G%s%1.0fD%s%1.0fP%s%1.0f%sT%s1',spaces(6),spaces(6),1,spaces(k),obj_num+1,spaces(p),num_lines+1,spaces(40),spaces(6));  % Nothing for now


% --- Create Spaces
function sp=spaces(n),

sp(1:n) = ' ';