function cnctoolbar(varargin),

global dcDebug

if nargin,
    feval(varargin{1}),
else
    loadtoolbar
end

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    if nargin
        fprintf(1,'  %s\n',varargin{1})
    end
end


% =====================================================
%                  Load CNC Toolbar
% =====================================================
function loadtoolbar

fig = gcf;
ax = get(fig,'currentaxes');
axis tight

ht = findall(fig,'type','uitoolbar');
if isempty(ht)
    ht = uitoolbar(fig);
else
%     rmbuttons(ht);
    delete(ht);
    ht = uitoolbar;
end
% set(ht,'tag','cnctoolbar')
color = 0.5*[1 1 1];
set(fig,'numbertitle','off','name','CNC Tool','color','k')
set(ax,'box','on','xcolor',color,'ycolor',color,'zcolor',color,'color','none')


% Set Default parameters
h = findobj(ax,'type','surface');
hs = get(h);
if length(hs) > 1,
    hs = hs(1);
end
% udata.htext= tex;
udata.curv = true;
udata.dir  = 'v';
udata.flip = true;
udata.D    = .5;
udata.step = .1;
udata.tol  = .1;
udata.X    = hs.XData;
udata.Y    = hs.YData;
udata.Z    = hs.ZData;
udata.surf = h;
udata.path = line(0,0,0);
udata.type = '';
xlim = [min(udata.X(:)) max(udata.X(:))];
xlim = xlim(1):diff(xlim)/4:xlim(2);
ylim = [min(udata.Y(:)) max(udata.Y(:))];
ylim = ylim(1):diff(ylim)/4:ylim(2);
zlim = [min(udata.Z(:)) max(udata.Z(:))];
set(ax,'xtick',xlim,'ytick',ylim,'ztick',zlim);

% Load icons
% load('icons.mat');
load('cncIcons.mat');
% (33) zoom_in.png
% (34) zoom_out.png
% (31) zigzag.png
% (27) spiral.png
% (14) flip.png
% (6) bullet_blue.png
% (13) eye.png
% (2) application_form.png
% (20) page_white_edit.png
% (9) cnc.png
% (11) color_wheel.png
%     'putdowntext(''rotate3d'',gcbo)'
%     'putdowntext(''pan'',gcbo)'
%     'putdowntext(''zoomout'',gcbo)'
%     'putdowntext(''zoomin'',gcbo)'

% Display Icons in Toolbar Menu
uitoggletool(ht,'CData',transp(icn(1).CData),'TooltipString','Zoom In','ClickedCallback',...
    'putdowntext(''zoomin'',gcbo)','tag','ZoomIn','separator','on')
uitoggletool(ht,'CData',transp(icn(2).CData),'TooltipString','Zoom Out','ClickedCallback',...
    'putdowntext(''zoomout'',gcbo)','tag','ZoomOut')

uitoggletool(ht,'CData',[],'TooltipString','Move','ClickedCallback',...
    'putdowntext(''zoomout'',gcbo)','tag','Move')
uitoggletool(ht,'CData',[],'TooltipString','Rotate','ClickedCallback',...
    'putdowntext(''rotate3d'',gcbo)','tag','Rotate')

uipushtool(ht,'CData',transp(icn(3).CData),'TooltipString','Milling Parameters','ClickedCallback',...
    'cnctoolbar(''setPush'')','tag','cnc.set','separator','on')
uipushtool(ht,'CData',transp(icn(4).CData),'TooltipString','Parallel Toolpath','ClickedCallback',...
    'cnctoolbar(''zigPush'')','tag','cnc.zigzag','separator','on')
uipushtool(ht,'CData',transp(icn(5).CData),'TooltipString','Spiral Toolpath','ClickedCallback',...
    'cnctoolbar(''spiralPush'')','tag','cnc.spiral')
uipushtool(ht,'CData',transp(icn(6).CData),'TooltipString','Flip Toolpath','ClickedCallback',...
    'cnctoolbar(''flipPush'')','tag','cnc.flip')
uipushtool(ht,'CData',transp(icn(7).CData),'TooltipString','Show Data Points on Toolpath','ClickedCallback',...
    'cnctoolbar(''markerPush'')','tag','cnc.marker')
uipushtool(ht,'CData',transp(icn(8).CData),'TooltipString','Show Curvature On/Off','ClickedCallback',...
    'cnctoolbar(''setCurv'')','tag','cnc.curve','separator','on')
uipushtool(ht,'CData',transp(icn(9).CData),'TooltipString','Verify Toolpath','ClickedCallback',...
    'cnctoolbar(''verifyPush'')','tag','cnc.verify')
uipushtool(ht,'CData',transp(icn(10).CData),'TooltipString','Post Toolpath to NC file','ClickedCallback',...
    'cnctoolbar(''postPush'')','tag','cnc.post','separator','on')

% uipushtool(ht,'CData',[],'TooltipString','Rotate Path Direction','ClickedCallback',...
%     'cnctoolbar(''dirPush'')')
% uipushtool(ht,'separator','on')
% tex=uicontrol('style','text','units','normalized','position',[.8 0,.19,.08],...
%     'horizontalalignment','right','foregroundcolor','w',...
%     'backgroundcolor',get(fig,'color'),'parent',fig);
% udata.htext = tex;
set(ax,'userdata',udata);

setPush

function icn = transp(icon)
if isa(icon,'uint8')
    icon = double(icon)/255;
end    
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

% =====================================================
%              Rotate Path direction
% =====================================================
function dirPush
udata = get(gca,'userdata');
if strcmp(udata.dir,'v')
    udata.dir = 'h';
else
    udata.dir = 'v';
end
set(gca,'userdata',udata)


% =====================================================
%                   Plot Gradient
% =====================================================
function curvPush


udata = get(gca,'userdata'); 
h = udata.surf;
x = udata.X;
y = udata.Y;
z = udata.Z;

if udata.curv
    warning off
    [g,m] = surfature(x,y,z);

    % m = 1./m;

    set(h,'CData',m,'facecolor','interp','edgealpha',.2);

    if udata.flip
        set(h,'CData',m,'edgealpha',.1);
    else
        set(h,'CData',-m,'edgealpha',.1);
    end

    set(gca,'Clim',[-1/udata.D,1/udata.D])
    % set(gca,'Clim',[-udata.D,udata.D])

    warning on
else
    set(h,'Cdata',z,'facecolor','interp','edgealpha',.1)
    set(gca,'Clim',[min(z(:)),max(z(:))])
end

function setCurv
ax = gca;
udata = get(ax,'userdata');
udata.curv = ~udata.curv;
set(ax,'userdata',udata)
curvPush

% =====================================================
%                   Post Toolpath
% =====================================================
function postPush
udata = get(gca,'userdata');
path = udata.path;
x = get(path,'xdata');
y = get(path,'ydata');
z = get(path,'zdata');

z = z - udata.D/2;

prompt={'XY Feedrate (in/min)','Z Feedrate (in/min)','Spindle Speed (RPM)'};
name='Milling Parameters';
numlines=1;
udata      = get(gca,'userdata');
defaultanswer={'30','5','1800'};
a=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(a)
    return
end

data.feedrate(1) = str2num(a{1});
data.feedrate(2) = str2num(a{2});
data.rpm         = str2num(a{3});
data.units = 'in';


[filename, pathname, filterindex] = uiputfile( ...
       {'*.nc',  'NC-files (*.nc)'}, ...
        'Save as');
if isempty(filename),
    return
end
filename = strrep(filename,'.nc','');
data.filename = filename;

fid = fopen([pathname filename '.nc'],'w');
if fid ~= -1,
    format = 'X%0.4f Y%0.4f Z%0.4f\n';
    cncpost(fid,format,x,y,z,data);
    if fclose(fid) ~= -1;
        pause(5)
        open([data.filename '.nc']);
    end
end


% =====================================================
%                   Flip Toolpath
% =====================================================
function flipPush
udata = get(gca,'userdata');
udata.flip = ~udata.flip;
set(gca,'userdata',udata);
curvPush
title('Warning: Toolpath needs to be regenerated','color','r')


% =====================================================
%               Set Toolpath Parameters
% =====================================================
function setPush,
prompt={'Tool Diameter','Stepover','Tolerance'};
name='Parameters';
numlines=1;
udata      = get(gca,'userdata');
% defaultanswer={'0.5','0.1','0.1'};
D = udata.D;
step = udata.step;
tol = udata.tol;
defaultanswer = {num2str(D),num2str(step),num2str(tol)};
a=inputdlg(prompt,name,numlines,defaultanswer);
if isempty(a)
    return
end

udata.D    = str2num(a{1});
udata.step = str2num(a{2});
udata.tol  = str2num(a{3});

if udata.D-D | udata.step-step | udata.tol-tol
    title('Warning: Toolpath needs to be regenerated','color','r')
end
ax = gca;
set(ax,'userdata',udata);
curvPush

% =====================================================
%                  Zigzag Toolpath
% =====================================================
function zigPush
udata = get(gca,'userdata');
udata.type = 'zigPush';
set(gca,'userdata',udata);
[x,y,z] = surf2path(udata.X,udata.Y,udata.Z,udata.D,udata.step,udata.tol,'zigzag',udata.flip,udata.dir);
% [x,y,z] = prepToolpath(x,y,z);
set(udata.path,'Xdata',x,'Ydata',y,'Zdata',z)
axis tight
title('');
% [L,n] = pathdetails(x,y,z)

% =====================================================
%                   Spiral Toolpath
% =====================================================
function spiralPush
udata = get(gca,'userdata');
udata.type = 'spiralPush';
set(gca,'userdata',udata);
[x,y,z] = surf2path(udata.X,udata.Y,udata.Z,udata.D,udata.step,udata.tol,'spiral',udata.flip);
% [x,y,z] = prepToolpath(x,y,z);
set(udata.path,'Xdata',x,'Ydata',y,'Zdata',z)
axis tight
title('');
% [L,n] = pathdetails(x,y,z)


% =====================================================
%                   Verify Toolpath
% =====================================================
function verifyPush
udata = get(gca,'userdata');
h = udata.path;
X = get(h,'XData'); 
Y = get(h,'YData'); 
Z = get(h,'ZData'); 

% X = X - min(X(:));
% Y = Y - min(Y(:));
% Z = Z - max(Z(:));
D = udata.D;

figure('numbertitle','off','name','Verify Toolpath','color','k');

set(gca,'visible','off');
% axis tight
% color = 0*[1 1 1];
% set(gca,'xgrid','off','ygrid','off','zgrid','off','box','on',...
%     'xcolor',color,'ycolor',color,'zcolor',color,'color','none');

% VERIFY TOOLPATH
h = verifyToolpath(X,Y,Z,D,500);

ax  = get(h,'parent');
axes(ax);
axis tight
color = .5*[1 1 1];
xlim = get(ax,'xlim'); ylim = get(ax,'ylim'); zlim = get(ax,'zlim');
set(ax,'xgrid','off','ygrid','off','zgrid','off','box','on',...
    'xcolor',color,'ycolor',color,'zcolor',color,'color','none',...
    'xtick',xlim,'ytick',ylim,'ztick',zlim,'visible','on');



function markerPush
ax = gca;
udata = get(ax,'userdata');
h = udata.path;
marker = get(h,'marker');
if strcmpi(marker,'none')
    set(h,'marker','.');
else
    set(h,'marker','none');
end


function rmbuttons(ht)
uit = findall(ht,'type','uitoggletool');
uip = findall(ht,'type','uipushtool');
a = ((strcmp(get(uit,'tag'),'Exploration.DataCursor')));
b = ((strcmp(get(uit,'tag'),'Annotation.InsertLegend')));
c = ((strcmp(get(uit,'tag'),'Annotation.InsertColorbar')));
d = ((strcmp(get(uit,'tag'),'Standard.EditPlot')));
delete(uit(a|b|c|d));
a = ((strcmp(get(uip,'tag'),'Standard.PrintFigure')));
b = ((strcmp(get(uip,'tag'),'Standard.SaveFigure')));
c = ((strcmp(get(uip,'tag'),'Standard.FileOpen')));
d = ((strcmp(get(uip,'tag'),'Standard.NewFigure')));
e = ((strcmp(get(uip,'tag'),'Plottools.PlottoolsOn')));
f = ((strcmp(get(uip,'tag'),'Plottools.PlottoolsOff')));
delete(uip(a|b|c|d|e|f));


function [L,n] = pathdetails(x,y,z)
ax = gca;
udata = get(ax,'userdata');
L = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
n = length(x);
set(udata.htext,'string',sprintf('Path Length: %1.0f\n# Points:   %1.0f',[L,n]))