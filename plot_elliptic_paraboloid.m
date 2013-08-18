function plot_elliptic_paraboloid(origin,ab,h,varargin)

p = inputParser;
p.addRequired('origin');
p.addRequired('ab');
p.addRequired('h');
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[1 0 0]);
p.addParamValue('opacity',0.2);
p.addParamValue('edgeopacity',0.2);
p.addParamValue('N',20);
p.parse(origin,ab,h,varargin{:})

O  = p.Results.origin;
ab = p.Results.ab;
h  = p.Results.h;
R  = p.Results.rotate;
N  = p.Results.N;
if length(N) == 1
  N = [N N];
end

if length(ab) == 1
  a = ab; b = ab;
else
  a = ab(1); b = ab(2);
end

nu = linspace(0,2*pi,N(1)); % row
u = linspace(0,h,N(2))';    % column
 
x = a*sqrt(u/h)*cos(nu);
y = b*sqrt(u/h)*sin(nu);
z = h-u*ones(1,N(1));

pos = rotation_matrix_zyx(R)*[x(:) y(:) z(:)]';
xx = reshape(pos(1,:),size(x,1),size(x,2));
yy = reshape(pos(2,:),size(y,1),size(y,2));
zz = reshape(pos(3,:),size(z,1),size(z,2));

surf(O(1)+xx,O(2)+yy,O(3)+zz,...
  'facecolor',p.Results.colour,...
  'facealpha',p.Results.opacity,...
  'edgealpha',p.Results.edgeopacity)