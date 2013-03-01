function plot_coord(origin,varargin)
%% plot_coord( origin , <opts> )
%
% Plots a 3D coordinate system origin.
% By default is aligned with the right-handed XYZ global coordinate system.
% Use optional arguments to customise:
%
% KEY    VALUE    DESCRIPTION
%
% 'length' L      Length of axes lines
% 'headratio' R      Ratio of head size to line length (0 < R < 1)
% 'rotate' [u v w] Rotate coordinates 'u' degrees around the X-axis
%                                     'v' degrees around the Y-axis
%                                     'w' degrees around the Z-axis
% 'labels' [true/false] Whether to print coordinate system labels
% 'index'  str    String index on labels (default str='1')
% 'linecolour' [R G B] Red-Green-Blue colour of axis lines
% 'headcolour' [R G B] Red-Green-Blue colour of arrowhead faces
% 'headopacity' C Opacity of arrowhead faces
%

%% Parse inputs:

p = inputParser;
p.addRequired('origin');
p.addOptional('rotate',[0 0 0]);
p.addOptional('length',0.05);
p.addOptional('headratio',1/3);
p.addOptional('arrowangle',25);
p.addOptional('index','1');
p.addOptional('labels',true);
p.addParamValue('linecolour',[1 1 1]);
p.addParamValue('headcolour',0.5*[1 1 1]);
p.addParamValue('headopacity',0.9);
p.parse(origin,varargin{:})

O      = p.Results.origin;
al     = p.Results.length;
aratio = p.Results.headratio;
ang    = p.Results.arrowangle;

r    = p.Results.rotate;
ni   = p.Results.index;
ecol = p.Results.linecolour;
col  = p.Results.headcolour;
opac = p.Results.headopacity;
labels_bool = p.Results.labels;

% Constant that should actually be optional input
nameshift = al/5;

%% Definition of a single axis (X)
%
% This is rotated to plot the other two in Y and Z.

ax = [al; 0; 0]; % axis end point
pxy = ax - aratio*al*[cosd(ang);  sind(ang); 0]; % one point of the arrowhead
pxz = ax - aratio*al*[cosd(ang); 0;  sind(ang)]; % one point of the other arrowhead
head1 = [ax pxy pxy.*[1; -1; 1]]; % arrowhead points in XY plane
head2 = [ax pxz pxz.*[1; 1; -1]]; % arrowhead points in XZ plane

% Rotation matrix (note inverse order of application)
R = Rz(r(3))*Ry(r(2))*Rx(r(1));

%% Plot

hold on

if labels_bool
  text(-nameshift+O(1),-nameshift+O(2),O(3),['O_',ni]);
end
plot_one_coord(O,R*ax,        R*head1,        R*head2,        ['x_',ni],[2*nameshift; 0; 0])
plot_one_coord(O,R*Rz(+90)*ax,R*Rz(+90)*head1,R*Rz(+90)*head2,['y_',ni],[0; nameshift; 0])
plot_one_coord(O,R*Ry(-90)*ax,R*Ry(-90)*head1,R*Ry(-90)*head2,['z_',ni],[0; 0; nameshift])

%% Nested functions

  function plot_one_coord(O,a,head1,head2,name,ns)
    plot_line(O,a)
    if labels_bool
      text(ns(1)+O(1)+a(1),ns(2)+O(2)+a(2),ns(3)+O(3)+a(3),name);
    end
    patcht(O,head1)
    patcht(O,head2)
  end

  function patcht(O,R)
    patch(O(1)+R(1,:),O(2)+R(2,:),O(3)+R(3,:),ecol,'facealpha',opac,'facecolor',col);
  end

end

%% Sub functions
%
% These could all also be anonymous functions.

function plot_line(O,r)
plot3(O(1)+[0 r(1)],O(2)+[0 r(2)],O(3)+[0 r(3)],'k');
end

% Rotation matrices:

function R = Rz(t)
R = [cosd(t) -sind(t) 0;
     sind(t)  cosd(t) 0;
     0        0       1];
end

function R = Ry(t)
R = [ cosd(t) 0  sind(t);
      0       1  0;
     -sind(t) 0  cosd(t)];
end

function R = Rx(t)
R = [1 0        0      ;
     0 cosd(t) -sind(t);
     0 sind(t)  cosd(t)];
end


% assert( all( Rx(90)*[0;1;0]==[ 0; 0; 1]) )
% assert( all( Rz(90)*[0;1;0]==[-1; 0; 0]) )
% assert( all( Ry(90)*[1;0;0]==[ 0; 0;-1]) )
% assert( all( Rz(90)*[1;0;0]==[ 0; 1; 0]) )
% assert( all( Rx(90)*[0;0;1]==[ 0;-1; 0]) )
% assert( all( Ry(90)*[0;0;1]==[ 1; 0; 0]) )