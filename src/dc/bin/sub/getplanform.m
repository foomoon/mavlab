function planform = getplanform(fig)

if ~nargin
    fig = gcf;
    mode = 'start';
end

if isstr(fig)
    mode = fig;
    fig = gcf;
end

switch mode
    case 'start'
        figure(fig);
        set(fig,'name','Select Leading Edge','numbertitle','off')
        title('Select Leading Edge');
        
        ax = get(fig,'currentaxes');
        hlines = findobj(ax,'type','line');
        set(hlines,'buttondownfcn','getplanform(''selectLedge'')')
        uiwait(fig)
    case 'selectLedge'
        h = gcbo;
        pform.Ledge.x = get(h,'xdata');
        pform.Ledge.x = pform.Ledge.x - min(pform.Ledge.x);
        maxX = max(pform.Ledge.x);
        pform.Ledge.x = pform.Ledge.x/maxX;
        pform.Ledge.y = get(h,'ydata')/maxX;
        set(fig,'userdata',pform,'name','Select Trailing Edge')
        title('Select Trailing Edge');
        
        ax = get(fig,'currentaxes');
        hlines = findobj(ax,'type','line');
        set(hlines,'buttondownfcn','getplanform(''selectTedge'')')
        set(h,'buttondownfcn','','linewidth',3,'color','r');
    case 'selectTedge'
        h = gcbo;
        set(h,'buttondownfcn','','linewidth',3,'color','r');
        pform = get(fig,'userdata');
        pform.Tedge.x = pform.Ledge.x;
        x = get(h,'xdata');
        y = get(h,'ydata');
        x = x-min(x);
        maxX = max(x);
        x = x/maxX;
        y = y/maxX;
        pform.Tedge.y = spline(x,y,pform.Ledge.x);
        set(fig,'userdata',pform,'name','Select Trailing Edge')
        uiresume
    otherwise
        error('Incorrect input string')
end

if nargout,
    planform = get(fig,'userdata');
end

% h1 = msgbox('Select leading Edge');
% 
% while ishandle(h1)
% end
% 
% h2 = msgbox('Select trailing Edge');
% 
% while ishandle(h2)
% end