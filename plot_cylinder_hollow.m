function plot_cylinder_hollow( O, r, R, h, t, varargin )
%%
p = inputParser;
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('N',20);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',0.5);
p.addParamValue('edgeopacity',0.5);
p.parse(varargin{:})

N = p.Results.N;
Rot = rotation_matrix_zyx(p.Results.rotate);
col   = p.Results.colour;
opac  = p.Results.opacity;
eopac = p.Results.edgeopacity;

theta1 = linspace(t(1),t(2),N);
theta2 = linspace(t(2),t(1),N);

p1x = [r*cosd(theta1) R*cosd(theta2)];
p1y = [r*sind(theta1) R*sind(theta2)];
p1z = zeros(1,2*N);
p2z = h*ones(1,2*N);

p1 = Rot*[p1x p1x(1);p1y p1y(1);p1z p1z(1)];
p2 = Rot*[p1x p1x(1);p1y p1y(1);p2z p2z(1)];

patch(O(1)+p1(1,:),O(2)+p1(2,:),O(3)+p1(3,:),col,'facealpha',opac,'edgealpha',eopac)
patch(O(1)+p2(1,:),O(2)+p2(2,:),O(3)+p2(3,:),col,'facealpha',opac,'edgealpha',eopac)

for ii = [1:N-1 N:2*N]
patch(...
  O(1)+[p1(1,ii) p2(1,ii) p2(1,ii+1) p1(1,ii+1)],...
  O(2)+[p1(2,ii) p2(2,ii) p2(2,ii+1) p1(2,ii+1)],...
  O(3)+[p1(3,ii) p2(3,ii) p2(3,ii+1) p1(3,ii+1)],col,'facealpha',opac,'edgealpha',eopac)  
end

end

