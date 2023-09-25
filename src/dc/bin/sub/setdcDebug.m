function setdcDebug(flag)
% SETDCDEBUG  Toggle debug mode for MAVLAB
% 
global dcDebug
try
    islogical(dcDebug);
    dcDebug = ~dcDebug;
catch    
    dcDebug = true;
end

if nargin
    switch lower(flag)
        case 'on'
            dcDebug = true;
        case 'off'
            dcDebug = false;
        otherwise
            error('setdcDebug: Input must either be "on" or "off"')
    end
end

if dcDebug
    fprintf(1,'  DC Debug [0N] OFF\n')
else
    fprintf(1,'  DC Debug  ON [0FF]\n')
end