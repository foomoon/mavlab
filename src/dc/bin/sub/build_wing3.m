function [varargout] = build_wing3(wing),
% BUILD_WING3 Creates surface data for wing
%   BUILD_WING3 constructs 3-D wing surface data from the parameters
%   defined in the struct-array, "wing".
%
%   [X,Y,Z]=build_wing3(wing)
%
%   Input
%       wing:    (struct) containing Wing data
%
%   Output
%       (X,Y,Z):    (2D-arrays) containg wing surface points
%
% See also: MAVLAB,QUICKWING,TAPER,SWEEP,DIHEDRAL,BLENDSHAPES

% Make Data Names Easier to Use (Because I'm Lazy)
span = wing.span;
chord = wing.chord;
theta = wing.dihedral;
polyhedral = wing.polyhedral;
sweep_angle = wing.sweep;
twist_angle = wing.twist;
curvy = wing.curvature;
curvetype = wing.curvetype;

% Don't ask why I stored the planform data this way
xl=wing.ledge(:,1);
yl=wing.ledge(:,2);
xt=wing.tedge(:,1);
yt=wing.tedge(:,2);

% Scale Normalized Wing Data to Correct Span & Chord
xl = xl*span/2;
xt = xt*span/2;
yl = yl*chord;
yt = yt*chord;

% Sweep Planform to Specified Angle
[xl,yl] = sweep(xl,yl,sweep_angle);
[xt,yt] = sweep(xt,yt,sweep_angle);

% Determine the Z-coordinates of planform


if size(wing.ledge,2) > 2 | curvetype == [],
    zcurve = curvy*(span)*wing.ledge(:,3)';
else
    zcurve = zoffset(xl,span/2,chord/2,curvy,curvetype);
end
zdihedral = dihedral(xl,theta,polyhedral);
ztotal = zdihedral+zcurve;
wing.ztotal = ztotal;

% Sort Airfoil Data
[x,foils] = blendfoils(wing);

if ~wing.edgeref
    % Shift Airfoils to top edge reference
    for i=1:length(wing.ledge),
        chord_length(:,i) = (wing.ledge(:,2)-wing.tedge(:,2));
    end
    chord_length = wing.chord*chord_length';
    fn = foils.*chord_length;
    maxz = max(max(fn));
    maxperfoil = max(fn);
    maxoffset = (maxz-maxperfoil);
    ztotal = ztotal + maxoffset;
end

rotation_axis = .25; % About 1/4 Chord

% [y_twist,z_twist] = twist(x,xl,xt,yl,yt,twist_angle,rotation_axis);

% Assemble Mesh
for j=1:length(x),          % chord direction
    for i=1:length(xl),     % span direction        
        % Calculate Spanwise Twist
        theta_span = twist_angle*abs( ( xl(i)/(xl(end)) )*(pi/180) );
        chord = abs(yl(i)-yt(i));
        twist_matrix = [cos(theta_span) sin(theta_span); -sin(theta_span) cos(theta_span)];
        
        y_foil(i,j) = -chord*(x(j,i)-rotation_axis); 
        z_foil(i,j) = foils(j,i)*chord;
        shift_forward(i,j) = yl(i) - rotation_axis*chord; 
        
        temp = twist_matrix*[y_foil(i,j); z_foil(i,j)] + [shift_forward(i,j); 0];
        y_twist(i,j) = temp(1);
        z_twist(i,j) = temp(2);
        
        X(i,j) = xl(i); % Don't Mess
        Y(i,j) = y_twist(i,j);
        Z(i,j) = ztotal(i) + z_twist(i,j); 
    end
end

% % Calculate Spanwise Twist
% % for i=1:length(xl),
% %     theta_span(i) = twist_angle*abs( ( xl(i)/(xl(end)) )*(pi/180) );
% %     chord(i) = abs(yl(i)-yt(i));
% %     twist_matrix(:,:,i) = [cos(theta_span(i)) sin(theta_span(i)); -sin(theta_span(i)) cos(theta_span(i))];
% % end
% % % Assemble Mesh
% % for j=1:length(x),          % chord direction
% %     for i=1:length(xl),     % span direction
% %         
% %         theta_span(i) = twist_angle*abs( ( xl(i)/(xl(end)) )*(pi/180) );
% %         chord(i) = abs(yl(i)-yt(i));
% %         twist_matrix(:,:,i) = [cos(theta_span(i)) sin(theta_span(i)); -sin(theta_span(i)) cos(theta_span(i))];
% %         
% %         y_foil(i,j) = -chord(i)*(x(j,i)-rotation_axis); 
% %         z_foil(i,j) = foils(j,i)*chord(i);
% %         shift_forward(i,j) = yl(i) - rotation_axis*chord(i); 
% %         
% %         temp = twist_matrix(:,:,i)*[y_foil(i,j); z_foil(i,j)] + [shift_forward(i,j); 0];
% %         y_twist(i,j) = temp(1);
% %         z_twist(i,j) = temp(2);
% %         
% %         X(i,j) = xl(i); % Don't Mess
% %         Y(i,j) = y_twist(i,j);
% %         Z(i,j) = ztotal(i) + z_twist(i,j); 
% %     end
% % end

% Get rid of NaN (Due to unselected portions of wing)
nanref = abs(isnan(Y)-1);
[i,j] = find(nanref(:,1));
X = X(i,:);
Y = Y(i,:);
Z = Z(i,:);

% Take care of units
if isfield(wing,'units'),
    switch wing.units
        case 'in',
            kscale = 1;
        case 'ft',
            kscale = 12;
        case 'cm',
            kscale = 0.393700787;
        case 'm',
            kscale = 39.3700787;
        case 'fur',
            kscale = 7920;
        otherwise
            kscale = 1;
    end
else
    kscale = 1;
end

% Scale Wing Based on Units
if nargout == 3,
    varargout{1} = X*kscale;
    varargout{2} = Y*kscale;
    varargout{3} = Z*kscale;
else
    out.X = X;
    out.Y = Y;
    out.Z = Z;
    varargout{1} = out;
end