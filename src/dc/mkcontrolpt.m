function out = mkcontrolpt(hg)

if isstr(hg)
    mode = hg;
else
    mode = 'start';
end

switch mode
    case 'start'
        n = getworkpln(hg,'alpha');  % Transparency
        x = getworkpln(hg,'position');
        color = getcolor(x)*(n) + [1 1 1]*(1-n);        
        h = line(x,0,0,'markerfacecolor',color,'marker','o',...
            'markersize',10,'clipping','off','markeredgecolor','k');
        bdown = 'mkcontrolpt(''click'')';
        set(h,'userdata',hg,'buttondownfcn',bdown);
        if nargout, out = h; end
    case 'click'
        fig = gcf;
        ud.down   = get(fig,'windowbuttondownfcn');
        ud.motion = get(fig,'windowbuttonmotionfcn');
        ud.up     = get(fig,'windowbuttonupfcn');
        ud.udata  = get(fig,'userdata');
        ud.hdot   = gcbo;
        bmotion   = 'mkcontrolpt(''move'')';        
        set(fig,'windowbuttonmotionfcn',bmotion,'userdata',ud);  
        try, set(ud.udata.hsurf,'visible','off'); end
    case 'move'
        fig = gcbo;
        ax  = get(fig,'currentaxes');
        ud  = get(fig,'userdata');
        pt  = get(ax,'currentpoint');
        pt  = pt(1);
        pt  = max(pt,0);
        pt  = min(pt,1);
        set(ud.hdot,'xdata',pt);
        hg = get(ud.hdot,'userdata');
        n = getworkpln(hg,'alpha');  % Transparency
        x = getworkpln(hg,'position');
        color = getcolor(x)*(n) + [1 1 1]*(1-n);
        mvworkpln(hg,pt);
        set(ud.hdot,'markerfacecolor',color);
        bup = 'mkcontrolpt(''up'')';
        set(fig,'windowbuttonupfcn',bup);
    case 'up'
        fig = gcbo;
        ud = get(fig,'userdata');
        try, set(ud.udata.hsurf,'visible','on'); end        
        set(fig,'userdata',ud.udata,'windowbuttondownfcn',ud.down,...
            'windowbuttonmotionfcn',ud.motion,'windowbuttonupfcn',ud.up); 
        try, workpln2srf; end
    otherwise
end



function color = getcolor(x)
col = colormap;
y = linspace(0,1,length(col));
[m,i]=min(abs(y-x));
color = col(i,:);

