function varargout = dcMain(varargin)
% DCMAIN M-file for dcMain.fig
%      DCMAIN, by itself, creates a new DCMAIN or raises the existing
%      singleton*.
%
%      H = DCMAIN returns the handle to a new DCMAIN or the handle to
%      the existing singleton*.
%
%      DCMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCMAIN.M with the given input arguments.
%
%      DCMAIN('Property','Value',...) creates a new DCMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcMain_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcMain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcMain

% Last Modified by GUIDE v2.5 18-Feb-2007 19:14:20


% BEGIN Custom code

if nargin & ishandle(varargin{1})
    h = varargin{1};
    x = varargin{2};
    scrollmotionfcn(h,x);
    return
end



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcMain_OpeningFcn, ...
                   'gui_OutputFcn',  @dcMain_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dcMain is made visible.
function dcMain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcMain (see VARARGIN)

% Choose default command line output for dcMain
handles.output = hObject;

% Create Toolbar
uitoolbar_ButtonDownFcn(handles);

if ~isempty(varargin) & isstruct(varargin{1})
    handles.wing = varargin{1};
else
    handles.wing = mkWing;
end


set(handles.figure1,'handlevisibility','on')

dcRotate(handles)

% Initialize Wing parameters into their proper figure fields
handles = setparams(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcMain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dcMain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on mouse press over axes background.
function planformAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to planformAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% msgbox('Planform Editor')
dcGetPlan

% --- Executes on mouse press over axes background.
function airfoilAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to airfoilAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% foil = dcGetFoil
% foil = airfoiltool
hsurf = findall(handles.mainAxes,'type','surface');
% wing = get(hsurf,'userdata');
hfoil = findall(hObject,'tag','hfoil');
dcGetFoil(hsurf,hfoil);

% --- Executes on mouse press over axes background.
function spanAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to spanAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcGetSpan

function icn = transp(icon)
c = squeeze(icon(1,1,:));
R  = icon(:,:,1);
G  = icon(:,:,2);
B  = icon(:,:,3);
RR = R == c(1);
GG = G == c(2);
BB = B == c(3);
RM = RR&GG&BB;
R(RM) = nan;
G(RM) = nan;
B(RM) = nan;
icn(:,:,1) = R;
icn(:,:,2) = G;
icn(:,:,3) = B;

function uitoolbar_ButtonDownFcn(hObject, eventdata, handles)
if nargin==1
    handles = hObject;
    tag = 'start';
else
    tag = get(hObject,'tag');
    hsurf = findall(handles.mainAxes,'type','surface');
    wing = get(hsurf,'userdata');
end

switch tag
    case 'start'
        ht = uitoolbar(handles.figure1);
        % Load icons
%         load('icons.mat');
        load('mavlabIcons')
        
        % Display Icons in Toolbar Menu
        uipushtool(ht,'CData',transp(icn(1).CData),'TooltipString','New Wing','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','newtool','separator','on')
        uipushtool(ht,'CData',transp(icn(2).CData),'TooltipString','Open Wing file','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','opentool')
        uipushtool(ht,'CData',transp(icn(3).CData),'TooltipString','Save Wing','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','savetool')
        uipushtool(ht,'CData',transp(icn(7).CData),'TooltipString','Export Wing','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','exporttool')
        uipushtool(ht,'CData',transp(icn(4).CData),'TooltipString','3D Wing View Window','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','showtool','separator','on')
        uipushtool(ht,'CData',transp(icn(5).CData),'TooltipString','CNC','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','milltool','separator','on')
        uipushtool(ht,'CData',transp(icn(6).CData),'TooltipString','AVL','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','avltool')
        uipushtool(ht,'CData',transp(icn(8).CData),'TooltipString','Documentation','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','doctool','separator','on')
        uipushtool(ht,'CData',transp(icn(9).CData),'TooltipString','Help','ClickedCallback',...
            'dcMain(''uitoolbar_ButtonDownFcn'',gcbo,[],guidata(gcbo))','tag','questtool')
        
    case 'avltool'
        dcAVL(wing);
    case 'questtool'
    case 'milltool'        
        dcCNC(wing)
    case 'savetool'
        [filename, pathname, filterindex] = uiputfile( ...
            {'*.mat','MAT-files (*.mat)'; ...
            '*.mdl','Models (*.mdl)'; ...
            '*.*',  'All Files (*.*)'}, ...
            'Save as', 'Untitled.mat');
        save([pathname filename],'wing');
    case 'opentool'
        [filename, pathname, filterindex] = uigetfile( ...
            {'*.mat','MAT-files (*.mat)'; ...
            '*.mdl','Models (*.mdl)'; ...
            '*.*',  'All Files (*.*)'}, ...
            'Save as', 'Untitled.mat');
        try
            f = load([pathname filename]);
            n = fieldnames(f);
            wing = getfield(f,n{1});
            closereq
            mavlab(wing)
        catch
            disp('Action Canceled')
        end
    case 'showtool'
        figure;
        dcSurf(wing);
    otherwise
end

function updateWing(hObject,handles)
hsurf = findall(handles.mainAxes,'type','surface');
wing = get(hsurf,'userdata');
wing.span = str2double(get(handles.span,'string'));
wing.chord = str2double(get(handles.chord,'string'));
wing.camber = str2double(get(handles.camber,'string'));
wing.tip = str2double(get(handles.zscale,'string'));
wing.twist = str2double(get(handles.twist,'string'));
wing.sweep = str2double(get(handles.sweep,'string'));
wing.dihedral = str2double(get(handles.dihedral,'string'));
wing.edgeref = ~get(handles.maxcambref,'value');
wing.mirror = get(handles.mirror,'value');
[X,Y,Z] = dcBuild(wing);
set(hsurf,'Xdata',X,'Ydata',Y,'Zdata',Z,'CData',Z)
set(hsurf,'userdata',wing)
axes(handles.mainAxes)
axis tight
guidata(hObject,handles);


function handles = setparams(handles)
wing = handles.wing;
set(handles.span,'string',num2str(wing.span));
set(handles.chord,'string',num2str(wing.chord));
set(handles.camber,'string',num2str(wing.camber));
set(handles.zscale,'string',num2str(wing.tip));
set(handles.twist,'string',num2str(wing.twist));
set(handles.sweep,'string',num2str(wing.sweep));
set(handles.dihedral,'string',num2str(wing.dihedral));
set(handles.maxcambref,'value',~wing.edgeref);
set(handles.mirror,'value',wing.mirror);
[X,Y,Z] = dcBuild(wing);
axes(handles.mainAxes)
% handles.surf = surface(wing.X,wing.Y,wing.Z);
% set(handles.surf,'parent',handles.mainAxes)
x = wing.foils.x;
y = wing.foils.y;
hfoil = line(x,y,'parent',handles.airfoilAxes,'tag','hfoil','clipping','off');
x = wing.planform.Ledge.x;
y = wing.planform.Ledge.z;
hspan = line(x,y,'parent',handles.spanAxes);
hspan(2) = line(-x,y,'parent',handles.spanAxes);
xl = wing.planform.Ledge.x;
yl = wing.planform.Ledge.y;
xt = wing.planform.Tedge.x;
yt = wing.planform.Tedge.y;
hplan(1) = line(xl,yl,'parent',handles.planformAxes,'tag','hledge');
hplan(2) = line(xt,yt,'parent',handles.planformAxes,'tag','htedge');
hplan(3) = line(-xl,yl,'parent',handles.planformAxes,'tag','hledgen');
hplan(4) = line(-xt,yt,'parent',handles.planformAxes,'tag','htedgen');
hplan(5) = line([-xl(end) -xt(end)],[yl(end) yt(end)],'parent',handles.planformAxes);
hplan(6) = line([xl(end) xt(end)],[yl(end) yt(end)],'parent',handles.planformAxes);

hh = [hfoil(:); hspan(:); hplan(:)];
set(hh,'linewidth',3,'color','g','hittest','off')
sideViews = [handles.planformAxes, handles.airfoilAxes, handles.spanAxes];
pos = get(sideViews,'position');

n = .1;
for i=1:3
   axes(sideViews(i)); axis tight
   xlim = get(sideViews(i),'xlim');
   ylim = get(sideViews(i),'ylim');
   xlim = xlim + diff(xlim)*[-n n];
   ylim = ylim + diff(ylim)*[-n n];
   set(sideViews(i),'xlim',xlim,'ylim',ylim)
end

handles.surf = surface(X,Y,Z,'parent',handles.mainAxes);
set(handles.surf,'userdata',handles.wing)
set(handles.surf,'facecolor','interp','edgealpha',.1)
set(handles.mainAxes,'visible','off')
light('parent',handles.mainAxes)

% Scroll Buttons
scr = [handles.spanscroll,handles.chordscroll,handles.camberscroll,...
          handles.tipscroll,handles.sweepscroll,handles.twistscroll,...
          handles.dihedralscroll];
hs= uiscrollbtn(scr,'dcMain');
load buttons
for i=1:length(scr)
    set(scr(i),'CData',buttons(i).img,'tooltip',buttons(i).tooltip);
end

function [flag] = dcRotateCallback(obj,event_obj)
objTag = get(obj,'Tag');
switch objTag
    case 'airfoilAxes'
        flag = true;
    case 'planformAxes'
        flag = true;
    case 'spanAxes'
        flag = true;
    otherwise
        flag = false;
end


function dcRotate(handles)
% Turn on rotate
h = rotate3d(handles.figure1);
set(h,'enable','on','Buttondownfilter',@dcRotateCallback)
% setAllowAxesRotate(h,handles.mainAxes,true)
setAllowAxesRotate(h,handles.airfoilAxes,false)
setAllowAxesRotate(h,handles.spanAxes,false)
setAllowAxesRotate(h,handles.planformAxes,false)


function scrollmotionfcn(hObject,val)
handles = guihandles(hObject);
fig = get(hObject,'parent');
tag = get(hObject,'tag');
% fprintf(1,'%s: ',tag)
% fprintf(1,'%1.6f\n',val)

switch tag
    case 'spanscroll'
        h = handles.span;        
        val = 1*val;
    case 'chordscroll'
        h = handles.chord;        
        val = 1*val;
    case 'camberscroll'
        h = handles.camber;        
        val = 1*val/100;
    case 'tipscroll'
        h = handles.zscale;        
        val = 1*val/10;
    case 'twistscroll'
        h = handles.twist;        
        val = 1*val;
    case 'sweepscroll'
        h = handles.sweep;        
        val = 1*val;
    case 'dihedralscroll'
        h = handles.dihedral;        
        val = 1*val;
    otherwise
end

v = str2double(get(h,'string'));
s = sprintf('%0.3f',v + val);
set(h,'string',s)


updateWing(hObject,handles)






