%% Hatze's elliptic octo-paraboloid
%
% I think.
% Not much of a head.

%% One quarter

a = 5;
b = 6;
c = 8;

N = 10;
yrange = linspace(0,b,N);
xmax = a*sqrt(1-(yrange/b).^2);

y = repmat(yrange,[N 1]);
x = nan(size(y));
for ii = 1:N
  x(ii,:) = linspace(xmax(ii),0,N)';
end

k = sqrt(1-(y/b).^2);
z = c*k.*(1-(x./a).^8);

figure(1); clf; hold on
surf(x,y,z)
view(3)

%% All eight "quadrants"

a = 5;
b = 5;
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