function plot_parabolic_plate(O,ab,h,varargin)

p = inputParser;
p.addRequired('O');
p.addRequired('ab');
p.addRequired('h');
p.addParamValue('N',15);
p.addParamValue('colour',[0 0.5 0.5]);
p.addParamValue('opacity',1);
p.addParamValue('edgeopacity',1);
p.parse(O,ab,h,varargin{:})

col = p.Results.colour;
opac = p.Results.opacity;
eopac = p.Results.edgeopacity;
n = p.Results.N;

yy = linspace(-ab(2),ab(2),n);

x = [1;1]*(O(1) + -ab(1)/ab(2)^2 * yy.^2 + ab(1));
y = [1;1]*(O(2) + yy);
z = O(3) + h * (0:1)' * ones(1,n);

hold on
surf(x,y,z,'edgealpha',eopac,'facecolor',col,'facealpha',opac)
patch(x(1,:),y(1,:),z(1,:),col,'facealpha',opac,'edgealpha',eopac)
patch(x(2,:),y(2,:),z(2,:),col,'facealpha',opac,'edgealpha',eopac)
patch(x(1,[1 1 1 1]),y(1,[1 end end 1]),z([1 1 2 2],1),col,'facealpha',opac,'edgealpha',eopac)