function plot_elliptic_plate(O,r,h,varargin)

p = inputParser;
p.addRequired('O');
p.addRequired('r');
p.addRequired('h');
p.addOptional('segment',[0 1]);
p.addParamValue('rotate',[0 0 0]);
p.addParamValue('N',25);
p.addParamValue('colour',[0.5 0 0.5]);
p.addParamValue('opacity',0.5);
p.addParamValue('edgeopacity',0.5);
p.parse(O,r,h,varargin{:})

R     = p.Results.rotate;
if numel(R) == 3
  R = rotation_matrix_zyx(R);
end

t     = p.Results.segment;
N     = p.Results.N;
col   = p.Results.colour;
opac  = p.Results.opacity;
eopac = p.Results.edgeopacity;

theta = linspace(t(1)*2*pi,t(2)*2*pi,N);

pp = R*[r(1)*cos(theta), r(1)*cos(theta);
        r(2)*sin(theta), r(2)*sin(theta);
        zeros(size(theta)), repmat(h,size(theta))];

PP = repmat(O,[1 2*N]) + pp;

shape.FaceColor = col;
shape.Vertices = PP';

shape.Faces = [1:N;N+(1:N)];
patch(shape,'facealpha',opac,'edgealpha',eopac)

shape.Faces = [1:N-1; 2:N; N+(2:N); N+(1:N-1)]';
patch(shape,'facealpha',opac,'edgealpha',eopac)
