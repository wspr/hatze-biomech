function person = segment_arm(person,S)

P = person.segment(S).origin + person.segment(S).offset(:);
R = person.segment(S).Rglobal;
N = person.segment(S).Ncalc;
PI = person.const.pi;

ind = 1:N;

L = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2;
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);
b1 = person.meas{S-1}.all(3)/2;  %shoulder measurement used in fortran code

person.segment(S).a = a;
person.segment(S).b = b;

gamma = person.resample(person,S,cellfun( @(x) x(person.sex), person.density.arm ));
gamma_0 = person.density.humerous(person.sex);

Q = P+person.segment(S).Rglobal*[0;0;-L];
person.segment(S+1).origin = Q;

%% Calculations

% volume
v      = PI*a.*b*L/N;
v_0    = 2/3*PI*((b1/2)^3);     %b1 or b(1)?  hatze '79 says b(1), fortran code says b1
volume = sum(v) + v_0;

% mass
m    = gamma.*v;
m_0  = gamma(1)*v_0;            %fortran code has gamma(1) instead of gamma_0 ?
mass = m_0 + sum(m);

% Mass centroid:
xc = 0;
yc = 0;
zcc = m_0*0.375*(b1/2) - sum(m.*(2*ind-1).*L)/20;
zc=zcc./mass;

% Moments of inertia:
I_x = m.*(3*(b.^2)+(L/N).^2)./12;
I_y = m.*(3*(a.^2)+(L/N).^2)./12;
I_z = m.*(a.^2+b.^2)./4;

% principal moments of inertia;
c=(2/5-9/64);
Ip_x = m_0*((c*((b1/2)^2))+((0.375*b1/2)-zc)^2) + sum(I_x) + sum(m.*(L*(2*ind-1)/20+zc).^2);
Ip_y = m_0*((c*((b1/2)^2))+((0.375*b1/2)-zc)^2) + sum(I_y) + sum(m.*(L*(2*ind-1)/20+zc).^2);
Ip_z = m_0*(2*(b1/2)^2)/5 + sum(I_z);

% principal moments of inertia w.r.t local systems origin
PIOX=Ip_x+mass*zc^2;
PIOY=Ip_y+mass*zc^2;
PIOZ=Ip_z;

%coordinates of origin of axes
OX= 0;
OY= 0;
%OZ=(at1+d_x)/cos(theta7); where at1, d_x and theta7 are calculations from
%the shoulder

person.segment(S).Rlocal   = person.segment(person.segment(S).prior).Rlocal';
person.segment(S).Rglobal  = (person.segment(S).Rlocal*person.segment(person.segment(S).prior).Rglobal);
person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot...

if person.plot || person.segment(S).plot

  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  for ii = ind
    ph = L-ii*L/N; % plate height
    plot_elliptic_plate(Q+R*[0;0;ph],[a(ii) b(ii)],L/N,opt{:},'rotate',R)
  end

  % the hemisphere
  plot_sphere(P, b1/2, 'longrange',[0 1],'rotate',R,opt{:})  %b1 or b(1)?

end

end
