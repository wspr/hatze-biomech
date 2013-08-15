function plot_elliptic_plate(O,r,h,varargin)

p = inputParser;
p.addRequired('O');
p.addRequired('r');
p.addRequired('h');
p.addOptional('segment',[0 1]);
p.addParamValue('N',25);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',1);
p.addParamValue('edgeopacity',1);
p.parse(O,r,h,varargin{:})

t     = p.Results.segment;
N     = p.Results.N;
col   = p.Results.colour;
opac  = p.Results.opacity;
eopac = p.Results.edgeopacity;

r1 = [r(1);r(1)];
r2 = [r(2);r(2)];
theta = linspace(t(1)*2*pi,t(2)*2*pi,N);

x = O(1) + r1*cos(theta);
y = O(2) + r2*sin(theta);
z = O(3) + h*(0:1)'*ones(1,N);

hold on
surf(x,y,z,'edgealpha',eopac,'facecolor',col,'facealpha',opac)
%patch(x(1,:),y(1,:),z(1,:),col,'facealpha',opac,'edgealpha',eopac)
%patch(x(2,:),y(2,:),z(2,:),col,'facealpha',opac,'edgealpha',eopac)