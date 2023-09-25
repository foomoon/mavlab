function varargout = dcGetSpan(varargin)
% DCGETSPAN M-file for dcGetSpan.fig
%      DCGETSPAN, by itself, creates a new DCGETSPAN or raises the existing
%      singleton*.
%
%      H = DCGETSPAN returns the handle to a new DCGETSPAN or the handle to
%      the existing singleton*.
%
%      DCGETSPAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCGETSPAN.M with the given input arguments.
%
%      DCGETSPAN('Property','Value',...) creates a new DCGETSPAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcGetSpan_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcGetSpan_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcGetSpan

% Last Modified by GUIDE v2.5 20-Feb-2007 07:50:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcGetSpan_OpeningFcn, ...
                   'gui_OutputFcn',  @dcGetSpan_OutputFcn, ...
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
% --- Executes just before dcGetSpan is made visible.
function dcGetSpan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcGetSpan (see VARARGIN)

% Choose default command line output for dcGetSpan
handles.output = [];

handles.Ledge = line([0 1],[0 0],'parent',handles.Laxes,'color','k');
% handles.Tedge = line([0 1],[0 0],'parent',handles.Taxes,'color','k');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcGetSpan wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dcGetSpan_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try
    for i=1:length(handles.output)
        varargout{1} = handles.output{i};
    end
catch
    Disp('No Data Output')
end
% disp('output function')
closereq

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
set(handles.Ledge,'Xdata',x,'Ydata',y,'parent',handles.Laxes)

% --- Executes on button press in lmodpush.
function lmodpush_Callback(hObject, eventdata, handles)
% hObject    handle to lmodpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gspline2(handles.Ledge);

% ------------------------------------------------------------------------
%%   Trailing Edge
% ------------------------------------------------------------------------
function tExpedit_Callback(hObject, eventdata, handles)
% hObject    handle to tExpedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tExpedit as text
%        str2double(get(hObject,'String')) returns contents of tExpedit as a double

% --- Executes on button press in tExppush.
function tExppush_Callback(hObject, eventdata, handles)
% hObject    handle to tExppush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = get(handles.tExpedit,'string');
x = linspace(0,1,20);
str = formatstr(str);
eval(str);
set(handles.Tedge,'Xdata',x,'Ydata',y,'parent',handles.Taxes)

% --- Executes on button press in tmodpush.
function tmodpush_Callback(hObject, eventdata, handles)
% hObject    handle to tmodpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gspline2(handles.Tedge);

% ------------------------------------------------------------------------
%%   Parameters
% ------------------------------------------------------------------------
% --- Executes on slider movement.
function chordgapSlider_Callback(hObject, eventdata, handles)
% hObject    handle to chordgapSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'value');
str = num2str(val);
set(handles.chordgapEdit,'string',str);



function chordgapEdit_Callback(hObject, eventdata, handles)
% hObject    handle to chordgapEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chordgapEdit as text
%        str2double(get(hObject,'String')) returns contents of chordgapEdit
%        as a double
val = str2double(get(hObject,'string'));
set(handles.chordgapSlider,'value',val);



% ------------------------------------------------------------------------
%%   OK/Cancel
% ------------------------------------------------------------------------
% --- Executes on button press in OKpush.
function OKpush_Callback(hObject, eventdata, handles)
% hObject    handle to OKpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gap = get(handles.chordgapSlider,'value');
xl = get(handles.Ledge,'Xdata');
yl = get(handles.Ledge,'Ydata');
% xt = get(handles.Tedge,'Xdata');
% yt = get(handles.Tedge,'Ydata');
[xl,yl]=processData(xl,yl);
handles.output{1} = [xl(:) yl(:)];
% handles.output{2} = [xt(:) yt(:)];
guidata(hObject,handles)
uiresume
% closereq

% --- Executes on button press in Cancelpush.
function Cancelpush_Callback(hObject, eventdata, handles)
% hObject    handle to Cancelpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume
% closereq


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

