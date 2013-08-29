function plot_trapzoidal_wedge(O,b,c,l,h,varargin)
% origin is the *centre* of the trapezoid.
% Straight edges are along X and perpendicular to Y. Height is Z.

if nargin == 0
  O = [0;0;0];
  b = 1;
  c = 2;
  l = 3;
  h = 0.5;
end

p = inputParser;
p.addRequired('O');
p.addRequired('a');
p.addRequired('b');
p.addRequired('c');
p.addRequired('h');
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',0.2);
p.addParamValue('edgeopacity',0.2);
p.parse(O,b,c,l,h,varargin{:})

vertex_l = [O O O O]'+[b/2 l/2 0; -b/2 l/2 0; -c/2 -l/2 0; c/2 -l/2 0];
vertex_h = [O O O O]'+[b/2 l/2 h; -b/2 l/2 h; -c/2 -l/2 h; c/2 -l/2 h];

trapz.Vertices = [vertex_l;vertex_h];
trapz.Faces = [1 2 3 4;
               5 6 7 8;
               1 2 6 5;
               2 3 7 6;
               1 4 8 5;
               3 4 8 7];

patch(trapz,...
  'facecolor',p.Results.colour,...
  'facealpha',p.Results.opacity,...
  'edgealpha',p.Results.edgeopacity)