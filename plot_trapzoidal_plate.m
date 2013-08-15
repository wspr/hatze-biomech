function plot_trapzoidal_plate(O,a,b,c,h)

% O  MUST  always be a 3x1 column vector

if nargin == 0
  O = [0;0;0];
  a = 1;
  b = 2;
  c = 3;
  h = 0.5;
end

vertex_l = [O O O O]'+[a/2 c/2 0; -a/2 c/2 0; -b/2 -c/2 0; b/2 -c/2 0];
vertex_h = [O O O O]'+[a/2 c/2 h; -a/2 c/2 h; -b/2 -c/2 h; b/2 -c/2 h];

trapz.Vertices = [vertex_l;vertex_h];
trapz.Faces = [1 2 3 4;
               5 6 7 8;
               1 2 6 5;
               2 3 7 6;
               1 4 8 5;
               3 4 8 7];

patch(trapz,'facecolor','red','facealpha',0.5)