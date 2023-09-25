function varargout = drawMesh(varargin)
% DRAWMESH M-file for drawMesh.fig
%      DRAWMESH, by itself, creates a new DRAWMESH or raises the existing
%      singleton*.
%
%      H = DRAWMESH returns the handle to a new DRAWMESH or the handle to
%      the existing singleton*.
%
%      DRAWMESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWMESH.M with the given input arguments.
%
%      DRAWMESH('Property','Value',...) creates a new DRAWMESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drawMesh_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drawMesh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help drawMesh

% Last Modified by GUIDE v2.5 27-Mar-2006 20:09:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawMesh_OpeningFcn, ...
                   'gui_OutputFcn',  @drawMesh_OutputFcn, ...
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


% --- Executes just before drawMesh is made visible.
function drawMesh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drawMesh (see VARARGIN)

% Choose default command line output for drawMesh
handles.output = hObject;

handles.wing = varargin{1};
handles.hbd = showPerimeter(handles.wing);
handles.dummy = figure('name','dummy','numbertitle','off','visible','off');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drawMesh wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drawMesh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in flexPush.
function flexPush_Callback(hObject, eventdata, handles)
newSpline
set(handles.BatPush,'enable','on')


% --- Executes on button press in BatPush.
function BatPush_Callback(hObject, eventdata, handles)
h = get(gca,'children');
hf = h(1);
if testForRigid(hf),
    wing.layup.boundary = bound;
    wing.layup.batton   = [];
    wing.layup.flex     = [];
    wing.layup.h0       = w;
    [cost,p,t,wing]     = meshWing(wing);
    return
end
closeRegion(hf);
% DRAW BATTON REINFORCED REGIONS
newLine(2);
handles.hf = hf;
guidata(hObject,handles);


% --- Executes on button press in OKPush.
function OKPush_Callback(hObject, eventdata, handles)
w = eval(get(handles.widthEdit,'string'));
hbd = handles.hbd;
hf = handles.hf;
wing = handles.wing;
%% DATA PROCESSING
% GET AXES CHILDREN (ALL LINES, PATCHES, ETC)
h = get(gca,'children');
% DELETE ANY SINGLE POINT LINES
h = rmSinglePoint(h);
% FIND 2 POINT LINES (BATTONS)
isBatton = find2Point(h);
% HANDLES FOR (BATTONS)
hb = h(isBatton);       
% WIDEN ALL BATTONS TO WIDTH (w)
wideBattons(hb,w);
% FIND ALL INTERSECTIONS OF CURVES IN THE AXIS
interSet(h);
% FIND CLOSEST POINT TO EACH INTERSECTION POINT OF BATTONS AND SWITCH THEM
% WITH THE INTERSECTION POINT
interSectBattons(hb);
% RESET AXIS
axis equal

%% MESH WING
wing.layup.boundary = formatBound(hbd);
wing.layup.batton   = formatBatton(hb);
wing.layup.flex     = formatFlex(hf,w);
wing.layup.h0       = w;
[cost,p,t,wing]     = meshWing(wing);

%% PLOT MESH
showMesh(p,t);


% --- Executes on selection change in MaterialsPop.
function MaterialsPop_Callback(hObject, eventdata, handles)
% hObject    handle to MaterialsPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns MaterialsPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MaterialsPop


function widthEdit_Callback(hObject, eventdata, handles)
% hObject    handle to widthEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthEdit as text
%        str2double(get(hObject,'String')) returns contents of widthEdit as a double



%% SUBROUTINES
function hbd = showPerimeter(wing);
%% BUILD WING (X/Y SWAP IS INTENTIONAL)
[X,Y,Z] = buildWing(wing);
%% DEFINE PERIMETER OF WING
bound1 = [Y(1,:)' X(1,:)'];
bound2 = flipud([Y(end,:)' X(end,:)']);
if all(bound1(end,:) == bound2(1,:)),
    bound2 = bound2(2:end,:);
end
bound = [bound1; bound2];
hbd = line(bound(:,1),bound(:,2),'linewidth',2,'color','k');
line(-bound(:,1),bound(:,2),'linewidth',2,'color','k','handlevisibility','off');


function isRigid = testForRigid(hf),
X = get(hf,'Xdata');
if length(X) == 1,
    isRigid = 1;
else
    isRigid = 0;
end

function B = formatBatton(hb),
n = length(hb);
B = [];
for i=1:n,
    X = get(hb(i),'Xdata');
    Y = get(hb(i),'Ydata');
    X = X(:);  Y = Y(:);
    B1 = [X(1:2)' Y(1:2)'];
    B2 = [X(3:4)' Y(3:4)'];
    B = [B; B1; B2];
end


function F = formatFlex(hf,w),
% GET LINE DATA
X = get(hf,'Xdata');
Y = get(hf,'Ydata');
% REMOVE REDUNDANT LAST POINT
X = X(1:end-1);
Y = Y(1:end-1);
% CLOSE THE SPLINE
[X,Y] = closedSpline(X,Y,w);
% ADJUST POINTS CLOSE TO INTERSECTIONS
udata = get(hf,'userdata');
xsect = udata.xsect;
[X,Y] = repClosePts(X,Y,xsect);
F = [X(:) Y(:)];


function F = formatBound(hb),
X = get(hb,'Xdata');
Y = get(hb,'Ydata');
udata = get(hb,'userdata');
xsect = udata.xsect;
[X,Y] = repClosePts(X,Y,xsect);
F = [X(:) Y(:)];


function [X,Y] = repClosePts(X,Y,xsect),
n = size(xsect,1);
for i=1:n,
   x = xsect(i,1);
   y = xsect(i,2);
   [m,ind] = min(sqrt((X-x).^2 + (Y-y).^2));
   X(ind) = x;
   Y(ind) = y;
end


function wideBattons(hb,w),
for i=1:length(hb),
   x = get(hb(i),'Xdata');
   y = get(hb(i),'Ydata');
   [X,Y] = wideLine(x,y,w);
   set(hb(i),'Xdata',X);
   set(hb(i),'Ydata',Y);   
end


function interSectBattons(hb),
for i=1:length(hb),
   X = get(hb(i),'Xdata');
   Y = get(hb(i),'Ydata');
   udata = get(hb(i),'userdata');   
   pt = udata.xsect;   x = pt(:,1);   y = pt(:,2);
   ref = [];
   for j=1:length(x),
       [m,ind] = min(sqrt((X-x(j)).^2 + (Y-y(j)).^2));
       ref(j,1) = ind;
   end
   X(ref) = x;
   Y(ref) = y;
   set(hb(i),'Xdata',X);
   set(hb(i),'Ydata',Y);   
end


function isBatton = find2Point(h),
X = get(h,'Xdata');
isBatton = zeros(size(h));
for i=1:length(h),
    if length(X{i}) == 2,
        isBatton(i) = 1;
    end
end
isBatton = find(isBatton);


function h = rmSinglePoint(h),
isSinglePt = zeros(size(h));
for i=1:length(h),
    xx = get(h(i),'xdata');
    if length(xx) <= 1, isSinglePt(i) = 1; end
end
isSinglePt = find(isSinglePt);
delete(h(isSinglePt));
% GET AXES CHILDREN AFTER DELETION OF SINGLE POINT LINES
h = get(gca,'children');


function closeRegion(hf),
x = get(hf,'Xdata');
y = get(hf,'Ydata');
set(hf,'Xdata',[x(:); x(1)])
set(hf,'Ydata',[y(:); y(1)]);





