%% Hatze's elliptic octo-paraboloid
%
% I think.
% Not much of a head.

%% All eight "quadrants"

a = 5;
b = 6;
c = 8;

kf = @(y)   sqrt(1-(y/b).^2);
zf = @(x,y) c*kf(y).*(1-(x./a).^8);

[xx,yy] = meshgrid(linspace(-a,a,40),linspace(-b,b,40));

zz = zf(xx,yy);

figure(1); clf; hold on
set(gcf,'color','white')

% offset by [1;2;3] -- change this to suit!!
surf(1+xx,2+yy,3+zz)
surf(1+xx,2+yy,3-zz)

axis equal
axis off
view(3)
zoom(3)

%% Maybe a better way?!?!

clear all
close all
clc

a = 5;
b = 6;
c = 8;

zf = @(x,y) c*sqrt(1-(y/b).^2).*(1-(x./a).^8);

N = 41;
t = linspace(-pi,pi,N);
d = linspace(0,1,N);

[tt,dd] = meshgrid(t,d);

x = dd*a.*cos(tt);
y = dd*b.*sin(tt);
z = zf(x,y);

figure(2); clf; hold on
set(gcf,'color','white')

opt = {'facecolor',[1 0 0],'facealpha',1};

surf(x,y,+z,opt{:})
surf(x,y,-z,opt{:})

x = a*sin(t);
y = b*cos(t);
z = zf(x,y);
surf([x;x],[y;y],[z;-z],opt{:})

axis equal
axis off
view(3)
zoom(2)