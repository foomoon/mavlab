function varargout = cnc(varargin)
% CNC M-file for cnc.fig
%      CNC, by itself, creates a new CNC or raises the existing
%      singleton*.
%
%      H = CNC returns the handle to a new CNC or the handle to
%      the existing singleton*.
%
%      CNC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CNC.M with the given input arguments.
%
%      CNC('Property','Value',...) creates a new CNC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cnc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cnc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: FLOWLINE, GUIDE, GUIDATA, GUIHANDLES, MAVLAB, MPSERVO1

% Edit the above text to modify the response to help cnc

% Last Modified by GUIDE v2.5 12-May-2006 22:34:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cnc_OpeningFcn, ...
                   'gui_OutputFcn',  @cnc_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cnc is made visible.
function cnc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cnc (see VARARGIN)

% try
% % winchangeicon(handles.figure1,which('MavLogo.ico'));
% name = get(hObject,'name');
% pause(.001)
% icon(101,name,fullfile(fileparts(which(mfilename)),'MavLogo.ico'));
% catch
% end

handles.wing=varargin{1};

% Choose default command line output for cnc
handles.output = hObject;

handles.figwidth = 1.5;

loadimg(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cnc wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function loadimg(handles),
set(handles.Background,'handlevisibility','on');
set(handles.Window,'handlevisibility','off');
img = imread('3.jpg');
imshow(img);
set(handles.Background,'handlevisibility','off');
set(handles.Window,'handlevisibility','on');


% --- Outputs from this function are returned to the command line.
function varargout = cnc_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function feedrate_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feedrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function feedrate_edit_Callback(hObject, eventdata, handles)
% hObject    handle to feedrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of feedrate_edit as text
%        str2double(get(hObject,'String')) returns contents of feedrate_edit as a double


% --- Executes on button press in post_pushbutton.
function post_pushbutton_Callback(hObject, eventdata, handles)
% contents = get(handles.popupmenu1,'String'); 
% tool = contents{get(handles.popupmenu1,'Value')};
% diameter = str2num(tool(1:3));
% scallop_height = str2num(get(handles.scallop_edit,'string'));
% handles.wing.border = str2num(get(handles.border_edit,'string'));
% if get(handles.vertical_radiobutton,'value'),
%     dir = 'v';
% elseif get(handles.circular_radiobutton,'value'),
%     dir = 'c';
% else
%     dir = 'h';
% end
% if findstr(lower(tool),'ball'),
%     tooltype = 'ball';
% else
%     tooltype = 'flat';
% end
% [X,Y,Z] = flowline(handles.wing,scallop_height,diameter,dir,tooltype);

[X,Y,Z,dia] = buildToolpath(handles);
feedrate = str2num(get(handles.feedrate_edit,'string'));
contents = get(handles.popupmenu1,'String'); 
tool = contents{get(handles.popupmenu1,'Value')};
[filename,pathname] = uiputfile('*.nc','Save NC File');
% try
    fullpath = [pathname filename];
    mpservo1(X,Y,Z,fullpath,tool,dia,feedrate);
% catch
%     disp('Error Creating Toolpath')
%     return
% end
close(handles.figure1);



% --- Executes on button press in verify_pushbutton.
function verify_pushbutton_Callback(hObject, eventdata, handles)
figure
[X,Y,Z,dia] = buildToolpath(handles);
ptsPerInch = 16;
L = lineLength(X(:),Y(:),Z(:));
numPts = round(L*ptsPerInch);
Path = curvspace([X(:),Y(:),Z(:)],numPts);
tic
verifyToolpath(Path(:,1),Path(:,2),Path(:,3),dia,[]);
toc

function [X,Y,Z,dia] = buildToolpath(handles),
scallop_height = str2num(get(handles.scallop_edit,'string'));
contents = get(handles.popupmenu1,'String'); 
tool = contents{get(handles.popupmenu1,'Value')};
dia = str2num(tool(1:4));
handles.wing.border = str2num(get(handles.border_edit,'string'));
if get(handles.vertical_radiobutton,'value'),
    dir = 'v';
elseif get(handles.circular_radiobutton,'value'),
    dir = 'c';
else
    dir = 'h';
end
if findstr(lower(tool),'ball'),
    tooltype = 'ball';
else
    tooltype = 'flat';
end
[X,Y,Z] = flowline(handles.wing,scallop_height,dia,dir,tooltype);


% --- Executes on button press in view_pushbutton.
function view_pushbutton_Callback(hObject, eventdata, handles)
[X,Y,Z,dia] = buildToolpath(handles);
pl = pathlength(X,Y,Z);
% h = figure;
% set(h,'Units','normalized',...
%     'Position',[.1,.1,.8,.8],...
%     'Color','k',...
%     'Name','Verify',...
%     'NumberTitle','off',...
%     'menubar','figure');
% try
% % winchangeicon(h,which('MavLogo.ico'));
% name = get(h,'name');
% pause(.001)
% icon(101,name,fullfile(fileparts(which(mfilename)),'MavLogo.ico'));
% catch
% end
line(X,Y,Z+dia/2,'color','c');
hold on

yyy = handles.wing.Y-max(max(handles.wing.Y));
yyy = yyy - min(min(yyy))/2;
surf(handles.wing.X,yyy,-(handles.wing.Z + max(max(-handles.wing.Z))),'edgecolor','none','FaceColor','interp');

line(X(1),Y(1),Z(1)+dia/2,'marker','o','color','g');
line(X(end),Y(end),Z(end)+dia/2,'marker','o','color','r');
[x,y,z]=mkbox(X,Y,Z,dia);
hold on
plot3(0,0,max(Z),'mo','linewidth',2)
% verifypath(X,Y,Z,dia/2);
hold off
axis equal
axis off
set(gca,'projection','perspective');
view([30,30]);


ceta = 2;  % Ellapsed Time Adjustment Coefficient
feedrate = str2num(get(handles.feedrate_edit,'string'));
stats = {['Estimated Run Time: ' ];...
         ['   ' sprintf('%1.0f',ceta*pl/feedrate) ' minutes'];...
         [' '];...
         ['Bounding Box: '];...
         ['   X: ' sprintf('%1.2f',x) ' in'];...
         ['   Y: ' sprintf('%1.2f',y) ' in'];...
         ['   Z: ' sprintf('%1.2f',z) ' in']};
% uicontrol('style','text',...
%     'string',stats,...    
%     'units','normalized',...    
%     'position',[.05,.06,.13,.18],...%     'backgroundcolor','none',...
%     'foregroundcolor','k',...
%     'horizontalalignment','left',...
%     'fontsize',8);
set(handles.specs,'string',stats);

rotate3D


% --- Executes during object creation, after setting all properties.
function scallop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scallop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function scallop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to scallop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scallop_edit as text
%        str2double(get(hObject,'String')) returns contents of scallop_edit as a double


function pl = pathlength(X,Y,Z),

% [s1,s2] = size(X);
% X = reshape(X,s1*s2,1);
% Y = reshape(Y,s1*s2,1);
% Z = reshape(Z,s1*s2,1);

dx = abs(diff(X));
dy = abs(diff(Y));
dz = abs(diff(Z));

x = sum(dx);
y = sum(dy);
z = sum(dz);

pl = sqrt(x^2 + y^2 + z^2);

function [x,y,z]=mkbox(X,Y,Z,dia),

r = dia/2;

lower_bot_left = [min(X)-r,min(Y)-r,min(Z)];
lower_bot_right = [max(X)+r,min(Y)-r,min(Z)];
lower_top_left = [min(X)-r,max(Y)+r,min(Z)];
lower_top_right = [max(X)+r,max(Y)+r,min(Z)];

upper_bot_left = [min(X)-r,min(Y)-r,max(Z)];
upper_bot_right = [max(X)+r,min(Y)-r,max(Z)];
upper_top_left = [min(X)-r,max(Y)+r,max(Z)];
upper_top_right = [max(X)+r,max(Y)+r,max(Z)];

lower_box = [lower_bot_left; lower_bot_right; lower_top_right; lower_top_left; lower_bot_left];
upper_box = [upper_bot_left; upper_bot_right; upper_top_right; upper_top_left; upper_bot_left];

leg1 = [lower_bot_left; upper_bot_left];
leg2 = [lower_bot_right; upper_bot_right];
leg3 = [lower_top_left; upper_top_left];
leg4 = [lower_top_right; upper_top_right];

hold on
line(lower_box(:,1),lower_box(:,2),lower_box(:,3),'color','w');
line(upper_box(:,1),upper_box(:,2),upper_box(:,3),'color','w');
line(leg1(:,1),leg1(:,2),leg1(:,3),'color','w');
line(leg2(:,1),leg2(:,2),leg2(:,3),'color','w');
line(leg3(:,1),leg3(:,2),leg3(:,3),'color','w');
line(leg4(:,1),leg4(:,2),leg4(:,3),'color','w');
hold off

x = abs(lower_bot_right-lower_bot_left);
y = abs(lower_top_left-lower_bot_left);
z = abs(upper_bot_left-lower_bot_left);
x = x(1);
y = y(2);
z = z(3);


% % --- Executes on button press in expand_pushbutton.
% function expand_pushbutton_Callback(hObject, eventdata, handles)
% % hObject    handle to expand_pushbutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% pos=get(handles.figure1,'position');
% pos(3) = handles.figwidth*pos(3);
% set(handles.figure1,'position',pos);
% handles.figwidth = handles.figwidth^-1;
% if handles.figwidth == 1.5,
%     set(hObject,'string','>>');
% else
%     set(hObject,'string','<<');
% end
% 
% guidata(hObject, handles);


% % --- Executes on button press in layout_pushbutton.
% function layout_pushbutton_Callback(hObject, eventdata, handles)
% % hObject    handle to layout_pushbutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% fig = figure('numbertitle','off','name','Layout Toolpath','color','w');
% try
% % winchangeicon(fig,which('MavLogo.ico'));
% name = get(fig,'name');
% pause(.001)
% icon(101,name,fullfile(fileparts(which(mfilename)),'MavLogo.ico'));
% catch
% end
% [X,Y,Z] = layoutpath(handles.wing);
% 
% yyy = handles.wing.Y-max(max(handles.wing.Y));
% yyy = yyy - min(min(yyy))/2;
% surf(handles.wing.X,yyy,-(handles.wing.Z ),'edgecolor','none','FaceColor','interp','facealpha',.7);
% light
% colormap([.4 .6 .8]);
% line(X,Y,-Z,'linewidth',2,'color','b');
% axis equal
% view(-30,30)
% 
% % uiwait(fig);
%  
% 
% feedrate = 30;
% diameter = .125;
% tool = '1/16 BALL ENDMILL';
% [filename,pathname] = uiputfile('*.nc','Save NC File');
% try
%     fullpath = [pathname filename];
%     mpservo1(X,Y,-Z,fullpath,tool,diameter,feedrate);
% catch
%     disp('Did Not Create Toolpath')
%     return
% end


% --- Executes on button press in vertical_radiobutton.
function vertical_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to vertical_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vertical_radiobutton
set(handles.horizontal_radiobutton,'value',0);
set(handles.circular_radiobutton,'value',0);

% --- Executes on button press in horizontal_radiobutton.
function horizontal_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to horizontal_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of horizontal_radiobutton
set(handles.vertical_radiobutton,'value',0);
set(handles.circular_radiobutton,'value',0);

% --- Executes on button press in circular_radiobutton.
function circular_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to circular_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of circular_radiobutton
set(handles.vertical_radiobutton,'value',0);
set(handles.horizontal_radiobutton,'value',0);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function border_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to border_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function border_edit_Callback(hObject, eventdata, handles)
% hObject    handle to border_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of border_edit as text
%        str2double(get(hObject,'String')) returns contents of border_edit as a double









% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


