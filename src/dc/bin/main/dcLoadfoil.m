function varargout = dcLoadfoil(varargin)
% DCLOADFOIL M-file for dcLoadfoil.fig
%      DCLOADFOIL, by itself, creates a new DCLOADFOIL or raises the existing
%      singleton*.
%
%      H = DCLOADFOIL returns the handle to a new DCLOADFOIL or the handle to
%      the existing singleton*.
%
%      DCLOADFOIL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DCLOADFOIL.M with the given input arguments.
%
%      DCLOADFOIL('Property','Value',...) creates a new DCLOADFOIL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dcLoadfoil_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dcLoadfoil_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dcLoadfoil

% Last Modified by GUIDE v2.5 04-Feb-2007 13:27:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dcLoadfoil_OpeningFcn, ...
                   'gui_OutputFcn',  @dcLoadfoil_OutputFcn, ...
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


% --- Executes just before dcLoadfoil is made visible.
function dcLoadfoil_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dcLoadfoil (see VARARGIN)

% Choose default command line output for dcLoadfoil
handles.output = hObject;

% ======================================================
% Load airfoil list into listbox
flistfile = which('dcflist.mat');
load(flistfile);
set(handles.foillistbox,'string',{flist.name});
handles.flist = flist;
% Set X-axis
h = line([0,1],[0,0],'color',[.6,.7,.7],'parent',handles.foilaxes);
% ======================================================

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dcLoadfoil wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dcLoadfoil_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in foillistbox.
function foillistbox_Callback(hObject, eventdata, handles)
% hObject    handle to foillistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot Airfoil
i = get(hObject,'Value');
x = handles.flist(i).x;
y = handles.flist(i).y;
h = line(x,y,'color','k','marker','.','parent',handles.foilaxes);

% Display Camber and postion
[camber,ind] = max(y);
camberpos = x(ind);
set(handles.camberedit,'string',num2str(camber*100))
set(handles.camberposedit,'string',num2str(camberpos*100))

% Display description
set(handles.descriptionedit,'string',handles.flist(i).description)

% Return to Non-modify mode
set(handles.camberedit,'enable','inactive');
set(handles.camberposedit,'enable','inactive');
set(handles.descriptionedit,'style','text');

% --- Executes during object creation, after setting all properties.
function foillistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foillistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function descriptionedit_Callback(hObject, eventdata, handles)
% hObject    handle to descriptionedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of descriptionedit as text
%        str2double(get(hObject,'String')) returns contents of descriptionedit as a double


% --- Executes during object creation, after setting all properties.
function descriptionedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to descriptionedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function camberedit_Callback(hObject, eventdata, handles)
% hObject    handle to camberedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camberedit as text
%        str2double(get(hObject,'String')) returns contents of camberedit as a double


% --- Executes during object creation, after setting all properties.
function camberedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camberedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modifypush.
function modifypush_Callback(hObject, eventdata, handles)
% hObject    handle to modifypush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.camberedit,'enable','on');
set(handles.camberposedit,'enable','on');
set(handles.descriptionedit,'style','edit');

% --- Executes on button press in OKpush.
function OKpush_Callback(hObject, eventdata, handles)
% hObject    handle to OKpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cancelpush.
function cancelpush_Callback(hObject, eventdata, handles)
% hObject    handle to cancelpush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function camberposedit_Callback(hObject, eventdata, handles)
% hObject    handle to camberposedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of camberposedit as text
%        str2double(get(hObject,'String')) returns contents of camberposedit as a double


% --- Executes during object creation, after setting all properties.
function camberposedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camberposedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Exportpushbutton.
function Exportpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Exportpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i = get(handles.foillistbox,'Value');
try
    x = handles.flist(i).x;
    y = handles.flist(i).y;
    f = [x(:) y(:)]';
catch
    return
end
[filename, pathname] = uiputfile( ...
    {'*.dat', 'Airfoil Format (*.dat)'; ...
    '*.*',                   'All Files (*.*)'}, ...
    'Save as');

if filename,
    filename = strrep(filename,'.dat','');

    fid = fopen([pathname filename '.dat'],'w');

    if fid~=-1,
        fprintf(fid,'%s\n',filename);
        fprintf(fid,'%f %f\n',f);
        fclose(fid);
    end
end

