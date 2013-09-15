function plot_hoof(O,a,b,h,varargin)
%% Ellipto-parabolic hoof

p = inputParser;
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('colour',[1 0 0]);
p.addParamValue('opacity',0.5);
p.addParamValue('edgeopacity',0.5);
p.addParamValue('N',20);
p.parse(varargin{:})

R = rotation_matrix_zyx(p.Results.rotate);

N = p.Results.N;
M = 2*N-1;

opt = {'facealpha',p.Results.opacity,'edgealpha',p.Results.edgeopacity};

theta = linspace(-pi,pi,M);
x3 = a*cos(theta);
y3 = b*sin(theta);
z3 = h-h/a^2*x3.^2;

pp = R*[x3, x3;
        y3, y3;
        zeros(1,M), z3];

PP = repmat(O,[1 2*M]) + pp;

shape.FaceColor = p.Results.colour;
shape.Vertices = PP';

% base
shape.Faces = 1:M;
patch(shape,opt{:})

% top
shape.Faces = M+[(1:M-1) ; (2:M) ; M-(2:M) ; M-(1:M-1)]';
patch(shape,opt{:})

% sides
shape.Faces = [1:M-1; 2:M; M+(2:M); M+(1:M-1)]';
patch(shape,opt{:})

end

