function plot_elliptic_paraboloid(origin,a,b,h,varargin)

p = inputParser;
p.addRequired('origin');
p.addRequired('a');
p.addRequired('b');
p.addRequired('h');
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[1 0 0]);
p.addParamValue('opacity',0.5);
p.addParamValue('N',20);
p.parse(origin,a,b,h,varargin{:})

O     = p.Results.origin;
a     = p.Results.a;
b     = p.Results.b;
h     = p.Results.h;

R = p.Results.rotate;
N = p.Results.N;

nu = linspace(0,2*pi,N); % row
u = linspace(0,h,N)';    % column
 
x = a*sqrt(u/h)*cos(nu);
y = b*sqrt(u/h)*sin(nu);
z = h-u*ones(1,N);

pos = rotation_matrix_zyx(R)*[x(:) y(:) z(:)]';
xx = reshape(pos(1,:),size(x,1),size(x,2));
yy = reshape(pos(2,:),size(y,1),size(y,2));
zz = reshape(pos(3,:),size(z,1),size(z,2));

surf(O(1)+xx,O(2)+yy,O(3)+zz)