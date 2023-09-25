function h = dcSurf(varargin)

global dcDebug

switch nargin,
    case 1
        wing = varargin{1};
        if ~isfield(wing,'X') | isempty(wing.X)
            wing = dcBuild(wing);
        end
        X = wing.X;
        Y = wing.Y;
        Z = wing.Z;
    case 3
        X = varargin{1};
        Y = varargin{2};
        Z = varargin{3};
    otherwise
        error('DCSURF: Incorrect # of inputs')
end

props = {'edgealpha',.2,'facecolor','b','facelighting','gouraud'};
hsurf = surface(X,Y,Z,props{:});
axis equal
axis tight
% set(gca,'box','on')
try set(gca,'ztick',[min(Z(:)) max(Z(:))]); catch end
try set(gca,'xtick',linspace(min(X(:)),max(X(:)),5)); catch end
% set(gcf,'color','w')
light

view(-130,30);

try 
    str = struct2str(specWing(wing));
    htxt = text('fontname','fixedwidth','string',str,'position',[0 -4 1]);
end


if nargout == 1,
    h = hsurf;
end

if dcDebug
    fprintf(1,'(%s)\n',mfilename)
end