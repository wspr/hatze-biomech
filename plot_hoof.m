function plot_hoof(origin,a,b,h,varargin)
%% Ellipto-parabolic hoof

p = inputParser;
p.addRequired('origin');
p.addRequired('a');
p.addRequired('b');
p.addRequired('h');
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[1 0 0]);
p.addParamValue('opacity',0.5);
p.addParamValue('N',20);
p.parse(origin,a,b,h,varargin{:})

O     = p.Results.origin;
a     = p.Results.a;
b     = p.Results.b;
h     = p.Results.h;

R = p.Results.rotate;
N = p.Results.N;
M = 2*N-1;

opt = {'facecolor',p.Results.colour,'facealpha',p.Results.opacity};

theta1 = linspace(-pi,pi,M);
x3 = a*cos(theta1);
y3 = b*sin(theta1);
z3 = h-h/a^2*x3.^2;

theta2 = linspace(-pi,0,N);
x4 = a*cos(theta2);
y4 = b*sin(theta2);
z4 = h-h/a^2*x4.^2;

surf(O(1)+[x3; x3],O(2)+[y3; y3],O(3)+[zeros(1,M); z3],opt{:}) % sides
%patch(O(1)+x3,O(2)+y3,O(3)+zeros(1,M),opt{:})    % base
surf(O(1)+[x4;x4],O(2)+[-y4; y4],O(3)+[z4;z4],[1 1 1],opt{:}) % top

end

