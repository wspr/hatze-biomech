function plot_elliptic_paraboloid(O,h,h2,a,b,n)

if nargin < 4
n = 20;
end

nu = linspace(0,2*pi,n); % row
u = linspace(0,h,n)';    % column
 
x = O(1)+a*sqrt(u/h)*cos(nu);
y = O(2)+b*sqrt(u/h)*sin(nu);
z = O(3)-h-h2+u*ones(1,n); %%h2 is the distance from origin of the paraboloid to the O12 or O15;

surf(x,y,z)