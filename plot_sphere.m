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
p.addParamValue('N',15);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',1);
p.addParamValue('edgeopacity',1);
p.parse(O,r,varargin{:})

col   = p.Results.colour;
opac  = p.Results.opacity;
eopac = p.Results.edgeopacity;
N     = p.Results.N;

%% Define geometry and coordinates

theta = linspace(-pi,pi,N);      % row
phi   = linspace(-pi/2,pi/2,N)'; % column

x1 = O(1) + r*cos(phi)*cos(theta);
y1 = O(2) + r*cos(phi)*sin(theta);
z1 = O(3) + r*sin(phi)*ones(1,N);

%% Plot

surf(x1,y1,z1,'edgealpha',eopac,'facecolor',col,'facealpha',opac)

