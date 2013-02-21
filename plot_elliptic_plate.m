function plot_elliptic_plate(O,r,h,varargin)


p = inputParser;
p.addRequired('O');
p.addRequired('r');
p.addRequired('h');
p.addParamValue('N',25);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',1);
p.addParamValue('edgeopacity',1);

p.parse(O,r,h,varargin{:})

col = p.Results.colour;
opac = p.Results.opacity;
n = p.Results.N;
eopac = p.Results.edgeopacity;

r1 = [r(1);r(1)];
r2 = [r(2);r(2)];
theta = (0:n)/n*2*pi;
sintheta = sin(theta);
sintheta(n+1) = 0;

x = O(1) + r1 * cos(theta);
y = O(2) + r2 * sintheta;
z = O(3) + h * (0:1)' * ones(1,n+1);


hold on
surf(x,y,z,'edgealpha',eopac,'facecolor',col,'facealpha',opac)
patch(x(1,:),y(1,:),z(1,:),col,'facealpha',opac,'edgealpha',eopac)
patch(x(2,:),y(2,:),z(2,:),col,'facealpha',opac,'edgealpha',eopac)