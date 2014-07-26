% plot_parabolic_wedge([0;0;0], [1 1], [2 2], 2.5,'humoffset',0.5,'humradius',0.3)

function plot_parabolic_wedge(O,ab1,ab2,h,varargin)

% if the first entry of FACE is TRUE, and HUMOFFSET and HUMRADIUS
% are defined, a hole is left in the wedge for the humerus.

p = inputParser;
p.addRequired('O');
p.addRequired('ab1');
p.addRequired('ab2');
p.addRequired('h');
p.addParamValue('skew',0);
p.addParamValue('drop',0);
p.addParamValue('face',[true true]);
p.addParamValue('humoffset',NaN);
p.addParamValue('humradius',NaN);
p.addParamValue('N',21);
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[0 0.5 0]);
p.addParamValue('opacity',0.4);
p.addParamValue('edgeopacity',1);
p.parse(O,ab1,ab2,h,varargin{:})

col = p.Results.colour;
opac = p.Results.opacity;
eopac = p.Results.edgeopacity;
n = p.Results.N;
t = p.Results.skew;
d = p.Results.drop;

y1 = linspace(-ab1(2),ab1(2),n);
y2 = linspace(-ab2(2),ab2(2),n);

x1 = d+(-ab1(1)/ab1(2)^2 * y1.^2 + ab1(1));
x2 = d+t+(-ab2(1)/ab2(2)^2 * y2.^2 + ab2(1));
z1 = zeros(1,n);
z2 = h*ones(1,n);

% hole
hole = true;
if any(isnan([p.Results.humoffset, p.Results.humradius]))
  hole = false;
end
pmin = [x1(1); y1(1); 0];
pmax = [x1(end); y1(end); 0];
pc = (pmin+pmax)/2;
c = pc + [p.Results.humoffset;0;0];
r = p.Results.humradius;

np = 20;
phi = -pi/2+linspace(0,2*pi,np);

cx1 = c(1) + r*sin(phi);
cy1 = c(2) + r*cos(phi);
cz1 = zeros(size(cx1));

shape = [pc, [x1;y1;z1], pc, [cx1; cy1; cz1]];

%% Transformation to local coordinates

R = rotation_matrix_zyx(p.Results.rotate);

p1 = R*[x1;y1;z1];
p2 = R*[x2;y2;z2];

x = O(1)+[p1(1,:);p2(1,:)];
y = O(2)+[p1(2,:);p2(2,:)];
z = O(3)+[p1(3,:);p2(3,:)];

% hole
tshape = R*shape;
xh = O(1)+tshape(1,:);
yh = O(2)+tshape(2,:);
zh = O(3)+tshape(3,:);

%% Plot everything

hold on

surf(x,y,z,'edgealpha',eopac,'facecolor',col,'facealpha',opac)

if p.Results.face(1)
  if hole
    patch(xh,yh,zh,col,'facealpha',opac,'edgealpha',eopac)
  else
    patch(x(1,:),y(1,:),z(1,:),col/2,'facealpha',opac,'edgealpha',eopac)
  end
end

if p.Results.face(2)
  patch(x(2,:),y(2,:),z(2,:),col,'facealpha',opac,'edgealpha',eopac)
end

patch(...
  [x(1,1) x(1,end) x(2,end) x(2,1)],...
  [y(1,1) y(1,end) y(2,end) y(2,1)],...
  [z(1,1) z(1,end) z(2,end) z(2,1)],...
  col,'facealpha',opac,'edgealpha',eopac)