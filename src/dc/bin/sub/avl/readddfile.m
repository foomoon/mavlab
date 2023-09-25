%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syntax: [dddata, control] = readddfile(filename)
% called from: avlab.m
%
% Translates AVL compatible derivative file into an AVLab compatible dddata struct
% and control struct array

function [dddata, control] = readddfile(filename)
TRUE = 1;
FALSE = 0;

% Reading the dddata struct: all fieldnames in the struct are lowercase due to matlab
% restrictions. When a variable would normally be represented as upper case the double
% lowercase letter has been used. The letter p represents the prime symbol:
%     i.e. clptot = cl'tot, clltot = CLtot

% define dddata struct
dddata = struct ( ...
    'alpha',  0, 'pb2v',  0, 'ppb2v', 0, ...
    'beta',   0, 'qc2v',  0, ...
    'mach',   0, 'rb2v',  0, 'rpb2v', 0, ...
    'cxtot',  0, 'cltot', 0, 'clptot', 0, ...
    'cytot',  0, 'cmtot', 0, ...
    'cztot',  0, 'cntot', 0, 'cnptot', 0, ...
    'clltot', 0, ...
    'cdtot',  0, ...
    'cdvis',  0, 'cdind', 0, ...
    'clff',   0, 'cdff',  0, ...
    'cyff',   0, 'e',     0, ...
    'clla',   0, 'cllb',  0, ...
    'cya',    0, 'cyb',   0, ...
    'cla',    0, 'clb',   0, ...
    'cma',    0, 'cmb',   0, ...
    'cna',    0, 'cnb',   0, ...
    'cllp',   0, 'cllq',  0, 'cllr', 0, ...
    'cyp',    0, 'cyq',   0, 'cyr', 0, ...
    'clp',    0, 'clq',   0, 'clr', 0, ...
    'cmp',    0, 'cmq',   0, 'cmr', 0, ...
    'cnp',    0, 'cnq',   0, 'cnr', 0 ...
);

fid = fopen(filename, 'r');

if fid == -1
    disp('Err: Cannot open file in readddfile')
    dddata = 0;
    control = 0;
    return
end

fgetl(fid);
% read first part of derivative file
while 1
    txt = fgetl(fid);
    if isempty(txt)
        continue
    end
    
    if txt == -1 % end of file
        disp ('Err: End of file reached before data fully collected in readddfile')
        fclose(fid);
        return
    end

    % this most likely means a blank line or is otherwise meaningless
    if length(txt) < 4
        continue
    end
    
    switch txt(3:7)
    case 'Alpha'
        [arg, txt] = getwordstr(txt); % filter name
        [arg, txt] = getwordstr(txt); % filter '='
        [arg, txt] = getwordstr(txt); dddata.alpha = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.pb2v = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.ppb2v = str2num(arg);
    case 'Beta '
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.beta = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.qc2v = str2num(arg);
    case 'Mach '
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.mach = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.rb2v = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.rbp2v = str2num(arg);
    case 'CXtot'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cxtot = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cltot = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clptot = str2num(arg);
    case 'CYtot'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cytot = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cmtot = str2num(arg);
    case 'CZtot'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cztot = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cntot = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cnptot = str2num(arg);
    case 'CLtot'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cltot = str2num(arg);
    case 'CDtot'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cdtot = str2num(arg);
    case 'CDvis'
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cdvis = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cdind = str2num(arg);
    case 'CLff '
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clff = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cdff = str2num(arg);
    case 'CYff '
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cyff = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.e = str2num(arg);
    case '-----' % file separator, move on to next section
        break;
    end % switch
end % while

% read second part of derivative file
while 1
    txt = fgetl(fid);
    if isempty(txt)
        continue
    end
    
    if txt == -1 % end of file
        disp ('Err: End of file reached before data fully collected in readddfile')
        fclose(fid)
        return
    end

    % this most likely means a blank line or is otherwise meaningless
    if length(txt) < 21
        continue
    end
    
    switch txt(16:21)
    case '   CLa'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clla = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cllb = str2num(arg);
    case '   CYa'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cya = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cyb = str2num(arg);
    case '   Cla'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cla = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clb = str2num(arg);
    case '   Cma'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cma = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cmb = str2num(arg);
    case '   Cna'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cna = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cnb = str2num(arg);
    case '   CLp'
        [arg, txt] = getwordstr(txt(15:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cllp = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cllq = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cllr = str2num(arg);
    case '   CYp'
        [arg, txt] = getwordstr(txt(16:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cyp = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cyq = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cyr = str2num(arg);
    case '   Clp'
        [arg, txt] = getwordstr(txt(16:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clp = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clq = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.clr = str2num(arg);
    case '   Cmp'
        [arg, txt] = getwordstr(txt(16:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cmp = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cmq = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cmr = str2num(arg);
    case '   Cnp'
        [arg, txt] = getwordstr(txt(16:end));
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cnp = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cnq = str2num(arg);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt);
        [arg, txt] = getwordstr(txt); dddata.cnr = str2num(arg);
        break; % This is the last instruction for the second part
    end % switch
end % while

template = struct ( ...
    'name' , 'null', ...
    'CLd'  , 0, ...
    'CYd'  , 0, ...
    'Cld'  , 0, ...
    'Cmd'  , 0, ...
    'Cnd'  , 0, ...
    'CDffd', 0, ...
    'ed'   , 0);

txt = fgetl(fid); % skip blank line
txt = fgetl(fid); % this line defines all the control surfaces
i = 0;
while 1
    [arg, txt] = getwordstr(txt); % control surface name
    name = arg;
    [arg, txt] = getwordstr(txt); % filter out control surface designation
    
    i = i + 1;
    ctrl = setfield(template, 'name', name);
    control(i) = ctrl;
    
    if size(txt) == 0
        break;
    end
end % while

% read third part of derivative file until EOF reached
while 1
    txt = fgetl(fid);
    if isempty(txt)
        continue
    end
    
    if txt == -1
        break; % end of file
    end

    % this most likely means a blank line or is otherwise meaningless
    if length(txt) < 19
        continue
    end

    switch txt(16:20)
    case '  CLd'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).CLd = str2num(arg);
        end % for
    case '  CYd'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).CYd = str2num(arg);
        end % for
    case '  Cld'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).Cld = str2num(arg);
        end % for
    case '  Cmd'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).Cmd = str2num(arg);
        end % for
    case '  Cnd'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).Cnd = str2num(arg);
        end % for
    case 'CDffd'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).CDffd = str2num(arg);
        end % for
    case '   ed'
        txt = txt(16:end);
        for j=1:1:i
            [arg, txt] = getwordstr(txt); % filter designation
            [arg, txt] = getwordstr(txt); % filter '='
            [arg, txt] = getwordstr(txt);
            control(j).ed = str2num(arg);
        end % for
    end
end % while

fclose(fid);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syntax: [arg, str] = getwordstr(str)
% called from: readrunfile, readddfile, readavlfile
%
% grabs the first argument from a string and returns the what is left over

function [arg, str] = getwordstr(str)

% return two blank strings if null-length input
if length(str) == 0
    str = '';
    arg = '';
    return
end

% strip leading spaces
i=1;
while isspace(str(i)) & i < length(str)
i=i+1;
end
str = str(i:length(str));

% capture target word
i=1;
while ~isspace(str(i)) & i < length(str)
    i=i+1;
end
arg = str(1:i);

if i < length(str)
    str = str(i:length(str));
else
    str = '';
end
return
