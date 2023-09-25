function installmavlab
% INSTALLMAVLAB  Type installmavlab to add mavlab to search path & toolbox
mavlabdir = uigetdir;

if ~mavlabdir
   return 
end

[pathstr,name] = fileparts(mavlabdir);

newpath = [matlabroot filesep 'toolbox' filesep name];
mkdir(newpath);
isok = copyfile(mavlabdir,newpath);

if isok    
    subpath = genpath(newpath);
    subpath = filterpathstr(subpath);
    for i=1:length(subpath)
        addpath(subpath{i},'-begin');
    end
    isok = ~savepath;
    if isok
        rehash toolboxcache
        disp('MAVLAB installed sucsessfully!')
    else
        warning('Problem saving path');
    end
else
    warning('Problem copying directory to toolpath')
end



function files = filterpathstr(str)
ind = findstr(str,';');
ind = ind(:);
ind0 = [1; ind(1:end-1) + 1];
ind =  ind - 1;

ind = [ind0 ind];

for i=1:size(ind,1)
    files{i} = str(ind(i,1):ind(i,2));
end

