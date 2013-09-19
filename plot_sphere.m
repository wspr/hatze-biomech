function plot_sphere(O,r,varargin)
%% plot_sphere( [x0 y0 z0], R, <opts> )
%
% - plots a sphere at [x0,y0,z0] of radius R.
% - optional arguments <opt> can be any combination of:
%
%   -------------------------------------------------------------
%   KEY            VALUE      DESCRIPTION
%   -------------------------------------------------------------
%   'N'            N          integer N number of segments to use
%   'colour'       [R G B]    Red-Green-Blue colour of each face
%   'opacity'      F          opacity of faces, F \in [0,1]
%   'edgeopacity'  F2         same for edges
%   -------------------------------------------------------------

%% Parse inputs

p = inputParser;
p.addRequired('O');
p.addRequired('r');
p.addParamValue('N',[15 5]);
p.addParamValue('latrange', [-1 1]);
p.addParamValue('longrange',[-1 1]);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',1);
p.addParamValue('edgeopacity',1);
p.addParamValue('rotate',[0 0 0]);
p.parse(O,r,varargin{:});

col   = p.Results.colour;
opac  = p.Results.opacity;
eopac = p.Results.edgeopacity;
N     = [p.Results.N(:); p.Results.N(:)]; % allow singleton input
rlong = p.Results.latrange;
rlat  = p.Results.longrange;
R = rotation_matrix_zyx(p.Results.rotate);

%% Define geometry and coordinates

theta = linspace(pi*rlong(1),pi*rlong(2),N(1));      % row
phi   = linspace(pi/2*rlat(1),pi/2*rlat(2),N(2))'; % column

x1 = r*cos(phi)*cos(theta);
y1 = r*cos(phi)*sin(theta);
z1 = r*sin(phi)*ones(1,N(1));

[f,v,c] = surf2patch(x1,y1,z1);

rv = repmat(O,[1 N(1)*N(2)]) + R*transpose(v);

%% Plot

patch(struct('Faces',f,'Vertices',transpose(rv),'FaceColor',col),...
  'edgealpha',eopac,'facealpha',opac);
