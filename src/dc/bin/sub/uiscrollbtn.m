function out = uiscrollbtn(varargin)

if nargin & isstr(varargin{1})
    mode = varargin{1};
    if isfunction(varargin{1})
        mode = 'start';
        fun = varargin{1};
    end    
elseif nargin & ishandle(varargin{1})
    % Recognize when input is a ui handle and convert object to
    % uiscrollbtn
    h = varargin{1};
    set(h,'style','pushbutton');
    mode = 'start';
    if nargin == 2
        fun = varargin{2};
    else
        fun = [];
    end
else
    mode = 'start'; 
    fun = [];
end






switch mode
    case 'start'
        if nargin > 2 & ~isempty(fun) & ~ishandle(varargin{1})
            h = uicontrol('style','pushbutton',varargin{2:end});
        elseif nargin & ~isempty(fun) & ~ishandle(varargin{1})
            h = uicontrol('style','pushbutton');
        elseif nargin & isempty(fun) & ~ishandle(varargin{1})
            h = uicontrol('style','pushbutton',varargin{:});
        elseif  ~ishandle(varargin{1})
            h = uicontrol('style','pushbutton');
        end
        
        set(h,'buttondownfcn','uiscrollbtn(''down'')','enable','inactive',...
            'userdata',fun)   
        if nargout
            out = h;
        end
    case 'down'       
        
        h = gcbo;
        scrnsize = get(0,'screensize');
        fig = getfig(h);
        hrot = rotate3d(fig);
        rotate = get(hrot,'enable');
        set(hrot,'enable','off');
        units = get(fig,'units');
        set(fig,'units','pixels');
        motion = get(fig,'windowbuttonmotionfcn');
        up = get(fig,'windowbuttonupfcn');
        udata = get(fig,'userdata');
        uu = get(h,'userdata');
        if isstruct(uu), 
            ud = uu;
        else
            ud.fun = uu; 
        end
        ud.hrot = hrot;
        ud.rotate = rotate;
        ud.scrnsize = scrnsize(3);
        ud.units = units;
        ud.motion = motion;
        ud.up = up;
        ud.udata = udata;
        ud.scroll = h;
        ud.firstpt = get(fig,'currentpoint');   
        ud.x = 0;
        set(fig,'userdata',ud);
        set(fig,'windowbuttonmotionfcn','uiscrollbtn(''motion'')',...
            'windowbuttonupfcn','uiscrollbtn(''up'')')
%         set(ud.scroll,'callback','')
    case 'motion'
        fig = gcbo;
        ud = get(fig,'userdata');
        firstpt = ud.firstpt;
        if isfield(ud,'lastpt')
            lastpt = ud.lastpt;
        else
            ud.lastpt = 0;
            lastpt = ud.lastpt;
        end
        ud.currpt = get(fig,'currentpoint');
        currpt = ud.currpt;
        dx = currpt(1) - firstpt(1);
        x = dx(1) + lastpt(1);
        ud.x = x;
        if ~isempty(ud.fun)
            feval(ud.fun,ud.scroll,x/ud.scrnsize)
        end
        set(fig,'userdata',ud);
        
    case 'up'
        fig = gcbo;
        ud = get(fig,'userdata');
        set(fig,'windowbuttonmotionfcn',ud.motion,'windowbuttonupfcn',...
            ud.up,'userdata',ud.udata,'units',ud.units)
        set(ud.hrot,'enable',ud.rotate);
        ud.lastpt = 0;%ud.x;
        set(ud.scroll,'buttondownfcn','uiscrollbtn(''down'')','userdata',ud)
    otherwise
end


function out = isfunction(str)
stats = functions(str2func(str));
if strcmp(stats.file,'')
    out = false;
else
    out = true;
end



function fig = getfig(h)

fig = get(h,'parent');
type = get(fig,'type');
if ~strcmpi(type,'figure')
    fig = getfig(fig);
end