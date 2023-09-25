function [name,coords]=scanfoil(filename,dlm),
% SCANFOIL  Airfoil Coordinate File Reader.
%   [name,coords]=SCANFOIL(filename,dlm) opens the file found in filename
%   and attempts to parse the file for airfoil coordinates delimited by the
%   character specified in dlm.  
%
%   SCANFOIL returns a name of the coordinates if specified in the first
%   lines, or returns 'default' otherwise.  It also returns the coordinates
%   in an m x 2 array, where m is the number of coordinates.
%   
%   See also SELECTSTATIONS, MAVLAB, FGETL


% Open file
fid=fopen(filename);

% Get airfoil name
name = fgetl(fid);
if ~isempty(str2num(name)),
    % No Name in file
    name = 'Default';
    fseek(fid,0,-1);
end

% Rename if first line is blank
if strcmp(name,''),
    name = 'Default';
end

% Determine Deliminator
if nargin < 2,
    [dlm,a] = finddlm(fid);
    fseek(fid,0,-1);
    if ~strcmp(name,'Default'),
        a=fgetl(fid);
    end
end

% Read in (X,Y) coordinates
i=1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    % if it is not a blank line
    if ~isempty(tline), 
        try
            coords(i,:)=sscanf(tline,['%f' dlm '%f'],[1,2]);
            i=i+1;
        catch
            [dlm,tline]=finddlm(fid);
            try
                coords(i,:)=sscanf(tline,['%f' dlm '%f'],[1,2]);
                i=i+1;
            catch
            end
        end
    end
end

% Close File
fclose(fid);




function [dlm,first_line]=finddlm(fid),

% Set Possible Deliminators
dlims = [',;'];
m = length(dlims);
dlm = ' ';
first_line = fgetl(fid);
% Search for first non empty line
while isempty(first_line),
    first_line = fgetl(fid);
end
n = length(first_line); 
% Search for deliminator
for i=1:n,
    for j=1:m,
        if strcmp(first_line(i),dlims(j)),
            dlm = dlims(j);
            break
        end
    end         
end    