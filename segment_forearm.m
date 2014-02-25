function person = segment_forearm(person,S)

P = person.segment(S).origin + person.segment(S).offset(:);
R = person.segment(S).Rglobal;
N = person.segment(S).Ncalc;
ind = 1:N;

L = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2;
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

gamma = person.resample(person,S,cellfun( @(x) x(person.sex), person.density.forearm ));

person.meas{S}.a = a;
person.meas{S}.b = b;

Q = P+person.segment(S).Rglobal*[0;0;-L];
person.segment(S+1).origin = Q;

%% Calculations 
% Mass
v = pi.*a.*b*L/N; 
m = gamma.*v; 

volume = sum(v);
mass = sum(m);

% Mass centroid (w.r.t original segment axes)
xc = 0;
yc = 0;
zc = -sum(m.*(ind-1/2)*L)/(N*mass);

% Moments of inertia
I_x = m.*(3*(b.^2)+(L/N).^2)/12;
I_y = m.*(3*(a.^2)+(L/N).^2)/12;
I_z = m.*(a.^2+b.^2)/4;

% principal moments of inertia
Ip_x = sum(I_x) + sum(m.*(L*(ind-1/2)/N+zc).^2);
Ip_y = sum(I_y) + sum(m.*(L*(ind-1/2)/N+zc).^2);
Ip_z = sum(I_z);

% principal moments of inertia w.r.t local systems origin
PIOX = Ip_x+mass*zc^2;
PIOY = Ip_y+mass*zc^2;
PIOZ = Ip_z;

%coordinates of origin of axes (origin of distal segment relative to local
%coordinate system of proximal segment)
OX=0;
OY=0;
OZ=-person.meas{S-1}.all(end);

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot

  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  for ii = ind
    ph = -ii*L/N; % plate height
    plot_elliptic_plate(P+R*[0;0;ph],[a(ii) b(ii)],L/N,opt{:},'rotate',R)
  end

end

end

