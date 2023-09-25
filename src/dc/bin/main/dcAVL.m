function varargout = dcAVL(varargin)
% DCAVL M-file for dcAVL.fig
%      DCAVL, by itself, creates a new DCAVL or raises the existing
%      singleton*.
%
%      H = DCAVL returns the handle to a new DCAVL or the handle to
%      the existing singleton*.
%
%      DCAVL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCAVL.M with the given input arguments.
%
%      DCAVL('Property','Value',...) creates a new DCAVL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcAVL_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcAVL_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcAVL

% Last Modified by GUIDE v2.5 15-Mar-2007 20:15:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcAVL_OpeningFcn, ...
                   'gui_OutputFcn',  @dcAVL_OutputFcn, ...
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


% --- Executes just before dcAVL is made visible.
function dcAVL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcAVL (see VARARGIN)

% Choose default command line output for dcAVL
handles.output = hObject;

if ~isempty(varargin) & isstruct(varargin{1})
    handles.wing = varargin{1};
end

loadLogo(handles);
handles = initParams(handles);
paramList_Callback(handles.paramList,[],handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcAVL wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dcAVL_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in paramList.
function paramList_Callback(hObject, eventdata, handles)
% hObject    handle to paramList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns paramList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from paramList
val = get(hObject,'value');
num = handles.paramsVal(val);
if length(val) > 1,
    set(handles.paramPanel,'title','Multiple Selection')
    set(handles.paramEdit,'string','')
else
    p = handles.paramsStr{val};
    set(handles.paramPanel,'title',p)    
    num = num2str(num);
    set(handles.paramEdit,'string',num)
end


% --- Executes during object creation, after setting all properties.
function paramList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function loadLogo(handles)

img = imread('avl.bmp');
image(img,'parent',handles.logoAxes);
set(handles.logoAxes,'xticklabel','','yticklabel','')

function handles = initParams(handles)
str = {'alpha','beta','pb/2V','qc/2V','rb/2V','CL','CDo','bank',...
        'elevation','heading','Mach','velocity','density','grav.acc',...
        'turn_rad.','load_fac','X_cg','Y_cg','Z_cg','mass','Ixx','Iyy',...
        'Izz','Ixy','Iyz','Izx'};

m = length(str);
handles.paramsStr = str;
p = zeros(m,1);


p(find(strcmp(str,'density'))) = 1.225;
p(find(strcmp(str,'grav.acc'))) = 9.81;

handles.paramsVal = p;

str = char(str);
str(:,end+1) = ' ';
str(:,end+1) = '=';
str(:,end+1) = ' ';
str(:,end+1) = ' ';
str(:,end+1) = ' ';


num = num2str(p,5);
str = [str num];

set(handles.paramList,'string',str);


function updateParams(hObject,handles)
val = get(handles.paramList,'value');

param = get(handles.paramEdit,'string');
param = strrep(param,' ','');
if isempty(param), param = 'nan'; end
param = eval(param);
handles.paramsVal(val) = param;
num = handles.paramsVal;
num = num2str(num);

str =  handles.paramsStr;
str = char(str);
str(:,end+1) = ' ';
str(:,end+1) = '=';
str(:,end+1) = ' ';
str(:,end+1) = ' ';
str(:,end+1) = ' ';

str = [str num];
set(handles.paramList,'string',str);
guidata(hObject,handles)

function paramEdit_Callback(hObject, eventdata, handles)
% hObject    handle to paramEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paramEdit as text
%        str2double(get(hObject,'String')) returns contents of paramEdit as a double
% val = get(handles.paramList,'value');
% str = get(hObject,'string');
% handles.paramsVal(val) = eval(str);
updateParams(hObject,handles);
% guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function paramEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OKpush.
function OKpush_Callback(hObject, eventdata, handles)
% hObject    handle to OKpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcDebug;
str = handles.paramsStr;
val = handles.paramsVal;
runData = struct;
for i=1:length(val)
    s = strrep(str{i},'/','3');
    s = strrep(s,'.','4');
    runData = setfield(runData,s,val(i));
end
if any(isnan(val))
    warndlg('Some incorrect values were entered','Warning')
    return
end
save runData runData
% Get Input params for avl
name = get(handles.nameEdit,'string');
wing = handles.wing;
Nspan = str2num(get(handles.NspanEdit,'string'));
Nchord = str2num(get(handles.NchordEdit,'string'));

% RUN AVL
[D] = avl(name,wing,runData,Nspan,Nchord);

% Format Output
str = struct2str(D);

% Display Output
figure('Numbertitle','off','name',sprintf('AVL Results: %s',name),...
       'toolbar','none','menubar','none');
uicontrol('style','listbox','string',str,'units','normalized','position',...
          [0 0 1 1],'fontname','fixedwidth');

% Display Debug info
if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    for i=1:size(str,1);
        fprintf(1,'%s\n',str(i,:));
    end
else
    clc;
end





function nameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to nameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameEdit as text
%        str2double(get(hObject,'String')) returns contents of nameEdit as a double


% --- Executes during object creation, after setting all properties.
function nameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NspanEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NspanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NspanEdit as text
%        str2double(get(hObject,'String')) returns contents of NspanEdit as a double


% --- Executes during object creation, after setting all properties.
function NspanEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NspanEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NchordEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NchordEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NchordEdit as text
%        str2double(get(hObject,'String')) returns contents of NchordEdit as a double


% --- Executes during object creation, after setting all properties.
function NchordEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NchordEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in plotPush.
function plotPush_Callback(hObject, eventdata, handles)
% hObject    handle to plotPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcDebug;
str = handles.paramsStr;
val = handles.paramsVal;
runData = struct;
for i=1:length(val)
    s = strrep(str{i},'/','3');
    s = strrep(s,'.','4');
    runData = setfield(runData,s,val(i));
end
if any(isnan(val))
    warndlg('Some incorrect values were entered','Warning')
    return
end

% Get Input params for avl
name = get(handles.nameEdit,'string');
wing = handles.wing;
Nspan = str2num(get(handles.NspanEdit,'string'));
Nchord = str2num(get(handles.NchordEdit,'string'));

alpha = linspace(0,20,10);

for i=1:length(alpha),
    runData.alpha = alpha(i);
    % RUN AVL
    [D] = avl(name,wing,runData,Nspan,Nchord);
    cl(i) = D.cltot;
    cm(i) = D.cmtot;
    cd(i) = D.cdtot;
    workbar(i/length(alpha),'AVL: Computing...');
end


figure('numbertitle','off','name',['Aerodynamics: ' name]);
subplot(2,2,1); plot(alpha,cl,'.-'); title('cl vs AOA'); 
subplot(2,2,2); plot(alpha,cm,'.-'); title('cm vs AOA'); 
subplot(2,2,3); plot(alpha,cd,'.-'); title('cd vs AOA'); 
subplot(2,2,4); plot(cd,cl,'.-'); title('cl vs cd'); 

out = [alpha(:) cl(:) cm(:) cd(:)]


