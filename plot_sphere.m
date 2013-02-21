function plot_sphere(O,r,varargin)
%%

p = inputParser;
p.addRequired('O');
p.addRequired('r');
p.addParamValue('N',15);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',1);

p.addParamValue('edgeopacity',1);
p.parse(O,r,varargin{:})

col = p.Results.colour;
opac = p.Results.opacity;
eopac = p.Results.edgeopacity;
n = p.Results.N;

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

x0 = cosphi*cos(theta);
y0 = cosphi*sintheta;
z0 = sin(phi)*ones(1,n+1);

x1 = x0*r+O(1);
y1 = y0*r+O(2);
z1 = z0*r+O(3);

surf(x1,y1,z1,'edgealpha',eopac,'facecolor',col,'facealpha',opac)

