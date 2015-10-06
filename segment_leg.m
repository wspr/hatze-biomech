function person = segment_leg(person,S)

P = person.segment(S).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;
N = person.segment(S).Ncalc;
ind = 1:N;

PI = person.const.pi;

%% Measurements

L = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2;
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

r = person.meas{S}.ankle/2;
r_a = 0.59*r;

gamma = person.resample(person,S,cellfun( @(x) x(person.sex), person.density.leg ));
gamma_b = person.density.ankle;


%% Calculations

% Volume
v      = PI*a.*b*L/N;
v_b    = 0.59*r^3*PI/2; % 0.59*r*PI/2 = 0.9268*r =~ 0.92*r as per Hatze and fortran code 
volume = sum(v) + 2*v_b;

% Mass
m    = gamma.*v;
m_b  = gamma_b*v_b; % ankle
mass = sum(m) + 2*m_b;

% Mass centroid:
xc = 0;
yc = 0;
zc = (-2*m_b*L - sum(m.*(2*ind-1)*L/20))./mass;

% Moments of inertia:
I_xi = m.*(3*b.^2+(L/N)^2)/12;
I_yi = m.*(3*a.^2+(L/N)^2)/12;
I_zi = m.*(a.^2+b.^2)/(4);  %fortran code has subroutine "ELZIN" which does:
% ... sum(gamma*a.*b.*u.^2)*L/251.3274 = m*u^2/(8*PI^2) = I_z/(2*PI) ...
%IF u^2/(2*PI^2)=a^2+b^2

% principal moments of inertia;
Ip_x = 2*m_b*(0.33*r^2+(L+zc)^2) + sum(I_xi + m.*(L*(2*ind-1)/20+zc).^2);
Ip_y = 2*m_b*(0.1859*r^2 + (L+zc)^2 + (a(N)+0.196*r)^2) + sum(I_yi + m.*(L*(2*ind-1)/20+zc).^2);
Ip_z = 2*m_b*(0.1859*r^2 + (a(N)+0.196*r)^2) + sum(I_zi);

% principal moments of inertia w.r.t local systems origin
PIOX=Ip_x+mass*zc^2;
PIOY=Ip_y+mass*zc^2;
PIOZ=Ip_z;

%coordinates of origin of axes
OX = 0;
OY = 0;
OZ = -person.meas{S-1}.length_long;


person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

Q = P+person.segment(S).Rglobal*[0;0;-L];
person.segment(S+1).origin = Q;

%% Plot

if person.plot || person.segment(S).plot

  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  for ii = ind
    ph = -ii*L/N; % plate height
    plot_elliptic_plate(P+R*[0;0;ph],[a(ii) b(ii)],L/N,opt{:},'rotate',R)
  end

  %% sideways paraboloids

  a = r;
  b = r;
  n = 10;

  nu = linspace(0,2*PI,n); % row
  u = linspace(0,r_a,n)';    % column

  x = a*sqrt(u/r_a)*cos(nu);
  y = b*sqrt(u/r_a)*sin(nu);
  z = u*ones(1,n);

  surf(Q(1)+z-a-r_a,Q(2)+y,Q(3)+x,'facealpha',person.segment(S).opacity(1),'edgealpha',person.segment(S).opacity(2),'facecolor',person.segment(S).colour)
  surf(Q(1)-z+a+r_a,Q(2)+y,Q(3)+x,'facealpha',person.segment(S).opacity(1),'edgealpha',person.segment(S).opacity(2),'facecolor',person.segment(S).colour)

end

end
