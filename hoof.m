%% Ellipto-parabolic hoof

% radii
a = 10;
b = 15;

% height
h = 20;

N = 40;

%%

figure(44); clf; hold on

theta = linspace(-pi,pi,N);

x = [a; a]*cos(theta);
y = [b; b]*sin(theta);
z = [h*ones(1,N); zeros(1,N)];
surf(x,y,z,'facecolor',[0 0 1],'facealpha',1,'edgealpha',0)

c = h/b^2;
xx = linspace(-a,a,N);
yy = linspace(-b,b,N);
zz = ones(N,1)*(h - c*yy.^2);
surf(xx,yy,zz,'facecolor',[1 0 0],'facealpha',1,'edgealpha',0)

x3 = x(1,:);
y3 = y(1,:);
z3 = h-h/a^2*x3.^2;
plot3(x3,y3,z3,'k.','markersize',50)

view(3)
camlight
axis equal

%%

figure(46); clf; hold on

opt = {'facecolor',[0 1 0],'facealpha',0.5};
surf(x,y,[zeros(1,N); z3],opt{:})
patch(x3,y3,zeros(1,N),opt{:})

theta = linspace(-pi,0,N);
x4 = a*cos(theta);
y4 = b*sin(theta);
z4 = h-h/a^2*x4.^2;

% fill towards "y=0" from both sides:
surf([x4;x4],[0*ones(1,N); y4],[z4;z4],[1 1 1],opt{:})
surf([x4;x4],[0*ones(1,N);-y4],[z4;z4],[1 1 1],opt{:})

view(3)
axis equal