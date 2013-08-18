function  plot_special_plate(origin,l,r,a,b,h,varargin)

p = inputParser;
p.addRequired('origin');
p.addRequired('l');
p.addRequired('r');
p.addRequired('a');
p.addRequired('b');
p.addRequired('h');
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[1 0 0]);
p.addParamValue('opacity',0.5);
p.addParamValue('edgeopacity',0.5);
p.addParamValue('N',9);
p.parse(origin,l,r,a,b,h,varargin{:})

O     = p.Results.origin;
l     = p.Results.l;
r     = p.Results.r;
a     = p.Results.a;
b     = p.Results.b;
h     = p.Results.h;
N     = p.Results.N;
R     = p.Results.rotate;

c1 = [ l/2-a; r-b; 0];
c2 = [-l/2+a; r-b; 0];

theta = linspace(0,pi/2,N);
p1 = repmat(c1,[1 N]) + [a*cos(theta);b*sin(theta);zeros(size(theta))];
theta = linspace(pi/2,pi,N);
p2 = repmat(c2,[1 N]) + [a*cos(theta);b*sin(theta);zeros(size(theta))];

bottom = [ [l/2;0;0] , p1 , p2 , [-l/2;0;0] ];
top = bottom;
top(3,:) = h;

bottom = rotation_matrix_zyx(R)*bottom + repmat(O,[1 size(bottom,2)]);
top = rotation_matrix_zyx(R)*top + repmat(O,[1 size(top,2)]);

opt = {'facecolor',p.Results.colour,'facealpha',p.Results.opacity,'edgealpha',p.Results.edgeopacity};

surf( [bottom(1,:);top(1,:)], [bottom(2,:);top(2,:)], [bottom(3,:);top(3,:)] , opt{:})
patch( bottom(1,:) , bottom(2,:) , bottom(3,:) , [0 0 0] , opt{:})
patch( top(1,:) , top(2,:) , top(3,:) , [0 0 0] , opt{:})
patch( [bottom(1,[1 end]) top(1,[end 1])] , ...
       [bottom(2,[1 end]) top(2,[end 1])] , ...
       [bottom(3,[1 end]) top(3,[end 1])] , [0 0 0] , opt{:})

end

