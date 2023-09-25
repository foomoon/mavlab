function varargout = dcGetFoil(varargin)
% DCGETFOIL M-file for dcGetFoil.fig
%      DCGETFOIL, by itself, creates a new DCGETFOIL or raises the existing
%      singleton*.
%
%      H = DCGETFOIL returns the handle to a new DCGETFOIL or the handle to
%      the existing singleton*.
%
%      DCGETFOIL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCGETFOIL.M with the given input arguments.
%
%      DCGETFOIL('Property','Value',...) creates a new DCGETFOIL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcGetFoil_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcGetFoil_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcGetFoil

% Last Modified by GUIDE v2.5 23-Apr-2007 15:33:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcGetFoil_OpeningFcn, ...
                   'gui_OutputFcn',  @dcGetFoil_OutputFcn, ...
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


% ------------------------------------------------------------------------
%%   Opening/Closing functions
% ------------------------------------------------------------------------
% --- Executes just before dcGetFoil is made visible.
function dcGetFoil_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcGetFoil (see VARARGIN)

if nargin > 3
    handles.hsurf = varargin{1};
    handles.hfoil = varargin{2};
    handles.wing = get(handles.hsurf,'userdata');
    x = handles.wing.foils.x(:,1);
    y = handles.wing.foils.y(:,1);
%     line(x,y,'color','r')
%     n = 11;
%     n = round((length(x)-2)/(n-2));    
%     x = [x(1); x(2:n:end-1); x(end)];
%     y = [y(1); y(2:n:end-1); y(end)];
else
    x = [0 1];
    y = [0 0];
end

% Choose default command line output for dcGetFoil
handles.output = [];

handles.Foil = line(x,y,'parent',handles.Laxes,'color','k','linewidth',2);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcGetFoil wait for user response (see UIRESUME)
% uiwait(handles.figure1);
gspline2(handles.Foil);

% --- Outputs from this function are returned to the command line.
function varargout = dcGetFoil_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    for i=1:length(handles.output)
        varargout{i} = handles.output{i};
    end
catch
    disp('No Data Output')
end
% disp('output function')
% figure(handles.figure1);
% closereq

% ------------------------------------------------------------------------
%%   Leading Edge
% ------------------------------------------------------------------------
function lExpedit_Callback(hObject, eventdata, handles)
% hObject    handle to lExpedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lExpedit as text
%        str2double(get(hObject,'String')) returns contents of lExpedit as a double


% --- Executes on button press in lExppush.
function lExppush_Callback(hObject, eventdata, handles)
% hObject    handle to lExppush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(handles.lExpedit,'string');
x = linspace(0,1,20);
str = formatstr(str);
eval(str);
set(handles.Foil,'Xdata',x,'Ydata',y,'parent',handles.Laxes)

% --- Executes on button press in lmodpush.
function lmodpush_Callback(hObject, eventdata, handles)
% hObject    handle to lmodpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gspline2(handles.Foil);
[filename, pathname] = uigetfile( ...
       {'*.dat;*.txt;*.pts', 'Airfoil Files (*.dat, *.txt, *.pts)'; ...
        '*.*',                   'All Files (*.*)'}, ...
        'Pick a file');
try
    [name,coords]=scanfoil([pathname filename]);
    x = coords(:,1);
    y = coords(:,2);
    set(handles.Foil,'Xdata',x,'Ydata',y,'ZData',x*0);
    gspline2(handles.Foil);
catch
    disp('Action Canceled');
end

% ------------------------------------------------------------------------
%%   Trailing Edge
% ------------------------------------------------------------------------
% function tExpedit_Callback(hObject, eventdata, handles)
% % hObject    handle to tExpedit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of tExpedit as text
% %        str2double(get(hObject,'String')) returns contents of tExpedit as a double

% % --- Executes on button press in tExppush.
% function tExppush_Callback(hObject, eventdata, handles)
% % hObject    handle to tExppush (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% str = get(handles.tExpedit,'string');
% x = linspace(0,1,20);
% str = formatstr(str);
% eval(str);
% set(handles.Tedge,'Xdata',x,'Ydata',y,'parent',handles.Taxes)

% % --- Executes on button press in tmodpush.
% function tmodpush_Callback(hObject, eventdata, handles)
% % hObject    handle to tmodpush (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% gspline2(handles.Tedge);

% ------------------------------------------------------------------------
%%   Parameters
% ------------------------------------------------------------------------
% --- Executes on slider movement.
% function chordgapSlider_Callback(hObject, eventdata, handles)
% % hObject    handle to chordgapSlider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% val = get(hObject,'value');
% str = num2str(val);
% set(handles.chordgapEdit,'string',str);


% 
% function chordgapEdit_Callback(hObject, eventdata, handles)
% % hObject    handle to chordgapEdit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of chordgapEdit as text
% %        str2double(get(hObject,'String')) returns contents of chordgapEdit
% %        as a double
% val = str2double(get(hObject,'string'));
% set(handles.chordgapSlider,'value',val);



% ------------------------------------------------------------------------
%%   OK/Cancel
% ------------------------------------------------------------------------
% --- Executes on button press in OKpush.
function OKpush_Callback(hObject, eventdata, handles)
% hObject    handle to OKpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gap = get(handles.chordgapSlider,'value');
x = get(handles.Foil,'Xdata');
y = get(handles.Foil,'Ydata');
% figure; plot(x,y); axis equal; 
% % xt = get(handles.Tedge,'Xdata');
% % yt = get(handles.Tedge,'Ydata');
% [xl,yl]=processData(xl,yl);
% handles.output{1} = [xl(:) yl(:)];
% % handles.output{2} = [xt(:) yt(:)];
% guidata(hObject,handles)
% uiresume(handles.figure1)
% % closereq
% pause(2)
% disp('this is a test')
hsurf = handles.hsurf;
hfoil = handles.hfoil;
wing = handles.wing;

p = curvspace([x(:) y(:)],51);
wing.foils.x = p(:,1);
wing.foils.y = p(:,2);
% wing.foils.x,return
% Z = Z*0;
[X,Y,Z] = dcBuild(wing);
set(hsurf,'Xdata',X,'YData',Y,'ZData',Z,'CData',Z,'userdata',wing);
set(hfoil,'Xdata',p(:,1),'Ydata',p(:,2));
% figure; surf(X,Y,Z,'edgecolor','none'); axis equal


% --- Executes on button press in Cancelpush.
function Cancelpush_Callback(hObject, eventdata, handles)
% hObject    handle to Cancelpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiresume
closereq


function str = formatstr(str)

str = strrep(str,'y','');
str = strrep(str,'=','');
str = strrep(str,'*','.*');
str = strrep(str,'/','./');
str = strrep(str,'^','.^');
str = strrep(str,'..','.');

str = ['y = ' str];



function [xl,yl,xt,yt]=processData(xl,yl,xt,yt,gap)

xl = xl/(max(xl)-min(xl));
% xt = xt/(max(xt)-min(xt));

yl = yl - yl(1);
% yt = yt - yt(1);
% yt = yt + yl(end) - yt(end) - gap;



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


