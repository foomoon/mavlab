function [data] = avl(name,wing,runData,Nspan,Nchord,varargin),

% Change to AVL directory
avldir = fileparts(which('avl26.exe'));
curdir = cd;
cd(avldir);

if ~isdir([cd filesep name])
    mkdir(name);
end
% Write AVL file from wing geometry
writeavlfile(name,wing,[Nspan Nchord]);
% Write AVL file from runData
writerunfile(name,runData);
% Write AVL input file
writeinpfile(name);
% Run AVL with the specified geometry
runavl(name);
% Read control derivative output
[data] = readddfile([name '\' name '_derivs.txt']);

% Change back to original directory
cd(curdir);




function runavl(name)
% Remove file extension if exists
name = strrep(name,'.run','');
name = strrep(name,'.avl','');
% Call avl from command prompt
system(['avl26 ' name '\' name ' < ' name '\' name '.inp']);


function writerunfile(name,runData)
% Remove file extension if exists
name = strrep(name,'.run','');
name = strrep(name,'.avl','');
% Open file for writing (overwrite if necessary)
fid = fopen([name '\' name '.run'],'w');

% Get AOA
AOA = runData.alpha;
% AOA = strrep(runData(1,:),'alpha     =','');
% AOA = strrep(AOA,' ','');
% AOA = str2double(AOA);

% Write run file data to file
fprintf(fid,'\n ---------------------------------------------\n');
fprintf(fid,' Run case  1:   %s (Created by MAVLAB)\n\n',name);
fprintf(fid,' alpha      ->  alpha       =   %f\n',AOA);
fprintf(fid,' beta       ->  beta        =   0.000000\n');
fprintf(fid,' pb/2V      ->  pb/2V       =   0\n');
fprintf(fid,' qc/2V      ->  qc/2V       =   0\n');
fprintf(fid,' rb/2V      ->  rb/2V       =   0\n\n');
% fprintf(fid,' aileron    ->  aileron     =   0.000000\n');
% fprintf(fid,' elevator   ->  elevator    =   0.000000\n');
% fprintf(fid,' rudder     ->  rudder      =   0.000000\n\n');

rd = struct2str(runData);
for i=1:length(rd)
    fprintf(fid,' %s\n',rd(i,:));
end
fprintf(fid,' visc CL_a = 0.000000\n');
fprintf(fid,' visc CL_u = 0.000000\n');
fprintf(fid,' visc CM_a = 0.000000\n');
fprintf(fid,' visc CM_u = 0.000000\n');
% Close file
fclose(fid);



function writeavlfile(name,wing,N)
% Remove file extension if exists
name = strrep(name,'.run','');
name = strrep(name,'.avl','');
% Open file for writing
fid = fopen([name '\' name '.avl'],'w');
% Do actual writing 
mav2avl(fid,wing,name,0,N);
% Close file
fclose(fid);


function writeinpfile(name)
% Remove file extension if exists
name = strrep(name,'.run','');
name = strrep(name,'.avl','');
% Open file for writing (overwrite if necessary)
fid = fopen([name '\' name '.inp'],'w');
% Command prompt inputs as they would be input in AVL
fprintf(fid,'oper\n');
fprintf(fid,'r\n');
fprintf(fid,'rm\n');
fprintf(fid,'0\n');
fprintf(fid,'x\n');
fprintf(fid,'o\n');
fprintf(fid,'p\n');
fprintf(fid,'F,F,F,T,F\n\n');

fprintf(fid,'w\n');
fprintf(fid,'%s\\%s_forces.txt\n',name,name);
fprintf(fid,'st \n');
fprintf(fid,'%s\\%s_derivs.txt\n',name,name);
fprintf(fid,'y\n\n');

fprintf(fid,'quit\n\n');

fprintf(fid,'quit\n');

% Close file
fclose(fid);


function str = struct2str(s)
% STRUCT2STR  Convert structure array to string
f = fieldnames(s);
n = length(f);

for i=1:n,
    temp = getfield(s,f{i});
    a = length(temp) > 1;
    b = isstruct(temp);
    c = iscell(temp);
    d = isempty(temp);
    if a | b | c | d
        temp = nan;
    end
%     temp
    val(i,1) = temp;
end

str1 = char(f);

for i=1:size(str1,1)
    str1(i,:) = strrep(str1(i,:),'3','/');
    str1(i,:) = strrep(str1(i,:),'4','.');
end

str1(:,end+1) = ' ';
str1(:,end+1) = '=';
str1(:,end+1) = ' ';
str1(:,end+1) = ' ';
str2 = num2str(val);
for i=1:size(str2,1)
    str2(i,:) = strrep(str2(i,:),'NaN','...');
end

str = [str1 str2];


%     'alpha',      6, ...       % angle of attack
%     'beta',       0, ...       % angle of sideslip
%     'pb2V',       0, ...       % nondimensionalised roll rate
%     'qc2V',       0, ...       % nondimensionalised pitch rate
%     'rb2V',       0, ...       % nondimensionalized yaw rate
%     'CL',         0, ...       % coefficient of lift
%     'CDo',        0, ...       % zero lift (form) drag coefficient
%     'bank',       0, ...       % bank (level = 0, right wing up is positive)
%     'elevation',  0, ...       % elevation via attitude (horizontal = 0)
%     'heading',    0, ...       % compass heading (East = 0)
%     'Mach',       0, ...       % mach number
%     'velocity',   15, ...       % ground velocity
%     'density',    1.225, ...% density of air
%     'gravacc',    9.81, ...    % acceleration due to gravity
%     'turn_rad',   0, ...       % turning radius
%     'load_fac',   1, ...       % load factor
%     'X_cg',       0.12, ...     % center of gravity ahead of geometry center
%     'Y_cg',       0, ...       % center of gravity right of geometry center
%     'Z_cg',       0, ...       % center of gravity below geometry center
%     'mass',       0.5, ...     % mass of aircraft
%     'Ixx',        0.0014, ...% moment of inertia xx-plane
%     'Iyy',        0.0079, ...% moment of inertia yy-plane
%     'Ixy',        0, ...       % moment of inertia xy-plane
%     'Iyz',        0, ...       % moment of inertia yz-plane
%     'Izx',        0, ...       % moment of inertia zx-plane
%     'viscCL_a',   0, ...       % _
%     'viscCL_u',   0, ...       % _
%     'viscCM_a',   0, ...       % _
%     'viscCM_u',   0);          % _