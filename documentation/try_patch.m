
clear all

NN = 15;
t = linspace(0,2*pi,NN);
x1 = cos(t);
x2 = cos(t);
y1 = sin(t);
y2 = sin(t);
z1 = zeros(size(x1));
z2 = ones(size(x1));

figure(77); clf; hold on
pp = rotation_matrix_zyx([30 0 0])*[x1, x2; y1, y2; z1, z2];

plot3(pp(1,:),pp(2,:),pp(3,:),'.','markersize',10)
for ii = 1:length(pp)
  text(pp(1,ii),pp(2,ii),pp(3,ii),['  ',num2str(ii)])
end
view(3)

shape.FaceColor = [1 0 0];
shape.Vertices = pp';

shape.Faces = [1:NN;NN+(1:NN)];
patch(shape,'facealpha',0.2)

shape.Faces = [1:NN-1; 2:NN; NN+(2:NN); NN+(1:NN-1)]';
patch(shape,'facealpha',0.2)

axis equal