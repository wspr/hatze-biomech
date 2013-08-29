function plot_parabolic_wedge(O,ab1,ab2,h,varargin)

p = inputParser;
p.addRequired('O');
p.addRequired('ab1');
p.addRequired('ab2');
p.addRequired('h');
p.addParamValue('t',0);
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
t = p.Results.t;

y1 = linspace(-ab1(2),ab1(2),n);
y2 = linspace(-ab2(2),ab2(2),n);

x1 = (-ab1(1)/ab1(2)^2 * y1.^2 + ab1(1));
x2 = t+(-ab2(1)/ab2(2)^2 * y2.^2 + ab2(1));
z1 = zeros(1,n);
z2 = h*ones(1,n);

if numel(p.Results.rotate) == 9
  R = p.Results.rotate;
elseif numel(R) == 3
  R = rotation_matrix_zyx(p.Results.rotate);
else
  error('Rotation input must be 3 cardan angles or the 3x3 rotation matrix itself.')
end
p1 = R*[x1;y1;z1];
p2 = R*[x2;y2;z2];

x = O(1)+[p1(1,:);p2(1,:)];
y = O(2)+[p1(2,:);p2(2,:)];
z = O(3)+[p1(3,:);p2(3,:)];

hold on
surf(x,y,z,'edgealpha',eopac,'facecolor',col,'facealpha',opac)
patch(x(1,:),y(1,:),z(1,:),col,'facealpha',opac,'edgealpha',eopac)
patch(x(2,:),y(2,:),z(2,:),col,'facealpha',opac,'edgealpha',eopac)
patch(...
  [x(1,1) x(1,end) x(2,end) x(2,1)],...
  [y(1,1) y(1,end) y(2,end) y(2,1)],...
  [z(1,1) z(1,end) z(2,end) z(2,1)],...
  col,'facealpha',opac,'edgealpha',eopac)