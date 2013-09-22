function plot_rect_prism(O,ab1,ab2,h,varargin)
% origin is the *centre* of the bottom face of the prism.
% Height is Z.


p = inputParser;
p.addRequired('O');
p.addRequired('ab1');
p.addRequired('ab2');
p.addParamValue('skew',[0 0]);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',0.2);
p.addParamValue('edgeopacity',0.2);
p.addParamValue('rotate',[0 0 0]);
p.parse(O,ab1,ab2,varargin{:})

R = rotation_matrix_zyx(p.Results.rotate);
s = p.Results.skew(1);
t = p.Results.skew(2);
a1 = ab1(1);
b1 = ab1(2);
a2 = ab2(1);
b2 = ab2(2);

vertex_l = [O O O O]+R*[a1/2 -a1/2 -a1/2 a1/2; b1/2 b1/2 -b1/2 -b1/2; 0 0 0 0];
vertex_h = [O O O O]+R*[s+a2/2 s-a2/2 s-a2/2 s+a2/2; t+b2/2 t+b2/2 t-b2/2 t-b2/2; h h h h];

trapz.Vertices = transpose([vertex_l,vertex_h]);
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