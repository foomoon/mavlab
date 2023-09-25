function val = getworkpln(h,param)

hpatch = findall(h,'tag','dcworkplane');
switch param
    case 'color'
        val = get(hpatch,'facecolor');
    case 'position'
        val = mean(get(hpatch,'xdata'));
    case 'alpha'
        val = get(hpatch,'facealpha');
    case 'nurb'
        hcurve = findall(h,'type','line');
        val = get(hcurve,'userdata');
    otherwise
end