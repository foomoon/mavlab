function [X,Y,Z]=prep4cnc(X,Y,Z)

global dcDebug

% Center Surface
Xshift = (-min(X(:)) - (max(X(:))-min(X(:)))/2);
Yshift = (-min(Y(:)) - (max(Y(:))-min(Y(:)))/2);
X = X + Xshift;
Y = Y + Yshift;

% Flip Surface
Z = -Z;

% Set max Z value to zero
Zshift = -max(Z(:));
Z = Z + Zshift;

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    fprintf(1,'  Surface Z Data inverted\n')
    fprintf(1,'  Surface shifted [%0.3f %0.3f %0.3f] \n',[Xshift Yshift Zshift])
end