function [xnew,ynew]=rmselfx(x,y)

%RMSELFX Self-intersections of a curve.
%
%    [XNEW,YNEW] = RMSELFX(X,Y) removes loops or kinks in a line in a fast 
%    and robust way.
%    The curve can be broken with NaNs or have vertical segments.
%    Segments of the curve involved in each of the self-interesections are
%    also provided.
%
%    Vectors X and Y are equal-length vectors of at least four points defining
%    the curve.
%
%    RMSELFX returns the x and y coordinates of the line resulting from the
%    removal of the loops/kinks.
%
%    This program uses the theory of operation of the file "Fast and Robust Curve
%    Intersections" submitted by Douglas M. Schwartz (intersections.m, F.Id: 11837).
%    It is a minor modification of Antoni J. Canos' "selfintersect.m"
%
%    Example of use
% 	 N=201;
% 	 th=linspace(-3*pi,4*pi,N);
% 	 R=1;
% 	 x=R*cos(th)+linspace(0,6,N);
% 	 y=R*sin(th)+linspace(0,1,N);
%    t0=clock;
%    [x0,y0]=rmselfx(x,y);
% 	 etime(clock,t0)
%    plot(x,y,'b',x0,y0,'.r');
% 	 axis ('equal'); grid

%
%    See also INTERSECTIONS.
%
%Version: 1.0, December 11, 2006
%Tested under MATLAB 6.5.0. R13.
%
% (c) Antoni J. Canos.
% ITACA. Techincal University of Valencia (Spain)
% Email:   ancama2@dcom.upv.es
%
% Modified: January 26, 2007
% University of Florida
% Daniel Claxton
% dclaxton@ufl.edu


% Input checks.
error(nargchk(2,2,nargin))
% x and y must be vectors with same number of points (at least 4 for self-intersection).
if sum(size(x) > 3) ~= 1 || sum(size(y) > 3) ~= 1 || ...
		length(x) ~= length(y)
	error('X and Y must be equal-length vectors of at least 4 points.')
end

x0=[];
y0=[];
segments=[];

% Two similar curves are firstly created.
x1=x; x2=x;
y1=y; y2=y;

x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;

dxy1 = diff([x1 y1]);
dxy2 = diff([x2 y2]);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
	repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
	repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
	repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
	repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
	repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Removing coincident and adjacent segments.
remove=find(abs(i-j)<2);
i(remove)=[];
j(remove)=[];

% Removing duplicate combinations of segments.
remove=[];
for ii=1:size(i,1)
	ind=find((i(ii)==j(ii:end))&(j(ii)==i(ii:end)));
	remove=[remove;ii-1+ind];
end
i(remove)=[];
j(remove)=[];

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  For some reason, the \ operation below is faster
% on my machine when A is sparse so we'll initialize a sparse matrix with
% the fixed values and then assign the changing values in the loop.
n = length(i);
T = zeros(4,n);
A = sparse([1 2 3 4],[3 3 4 4],-1,4,4,8);
B = -[x1(i) x2(j) y1(i) y2(j)].';
index_dxy1 = [1 3];  %  A(1) = A(1,1), A(3) = A(3,1)
index_dxy2 = [6 8];  %  A(6) = A(2,2), A(8) = A(4,2)

% Loop through possibilities.  Set warning not to trigger for anomalous
% results (i.e., when A is singular).
warning_state = warning('off','MATLAB:singularMatrix');
try
	for k = 1:n
		A(index_dxy1) = dxy1(i(k),:);
		A(index_dxy2) = dxy2(j(k),:);
		T(:,k) = A\B(:,k);
	end
	warning(warning_state)
catch
	warning(warning_state)
	rethrow(lasterror)
end

% Find where t1 and t2 are between 0 and 1 and return the corresponding x0
% and y0 values.  Anomalous segment pairs can be segment pairs that are
% colinear (overlap) or the result of segments that are degenerate (end
% points the same).  The algorithm will return an intersection point that
% is at the center of the overlapping region.  Because of the finite
% precision of floating point arithmetic it is difficult to predict when
% two line segments will be considered to overlap exactly or even intersect
% at an end point.  For this algorithm, an anomaly is detected when any
% element of the solution (a single column of T) is a NaN.

in_range = T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1;
anomalous = any(isnan(T));
if any(anomalous)
	ia = i(anomalous);
	ja = j(anomalous);
	% set x0 and y0 to middle of overlapping region.
	T(3,anomalous) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
		min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1))))/2;
	T(4,anomalous) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
		min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1))))/2;
	x0 = T(3,in_range | anomalous).';
	y0 = T(4,in_range | anomalous).';
	i=i(in_range | anomalous);
	j=j(in_range | anomalous);
else
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	i=i(in_range);
	j=j(in_range);
end

segments=sort([i,j],2);

s = segments;
xnew = x;
ynew = y;

if size(s,1) < 1,
    return
end

if diff(s(1,:)) < 3,
    %     return
    s(1,:) = [];
    x0(1) = [];
    y0(1) = [];
end

if size(s,1) > 1 & s(1,2) > s(2,1), 

% if any(s(1:end-1,2)>s(2:end,1))
    n = length(x);
    
%     fst = find(diff(s,1,2) > length(x)/2);
%     if length(fst)>1, [m,i] = max(diff(s(fst,:),1,2)); fst = fst(i); end
%     sfst = s(fst,:)
%     s(1:fst,:) = [];
%     x0(1:fst) = [];
%     y0(1:fst) = [];
%     s = [[0 sfst(1)]; s; [sfst(2) n]];
%     x0 = [x0(fst); x0; x0(fst)];
%     y0 = [y0(fst); y0; y0(fst)];
    
    first = s(1,:); 
    s(1,:) = [0 first(1)];
    s(end+1,:) = [first(2) n];
    x0 = [x0; x0(1)];
    y0 = [y0; y0(1)];
end


seg = []; 
for i=1:size(s,1), 
    seg = [seg s(i,1)+2:s(i,2)]; 
end
xnew(s(:,1)+1) = x0;
ynew(s(:,1)+1) = y0;
xnew(seg)=[]; ynew(seg) = [];