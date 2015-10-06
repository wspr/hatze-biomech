function person = segment_head_neck(person,S)

%% Head_neck

P = person.segment(S).origin+person.segment(S).offset;
R = person.segment(S).Rglobal;

i_m = person.sex;

head_width  = person.meas{S}.all(1);
head_depth  = person.meas{S}.all(2);
head_height = person.meas{S}.all(3);
neck_height = person.meas{S}.all(4);

l = person.meas{S-1}.length;
PI = person.const.pi;

%% Densities:

gamma_e = person.density.head(i_m);
gamma_c = person.density.neck(i_m);


%% Measurements

a = head_width/2;
b = 1.08*head_depth/2;
c = 1.04*head_height/2;

a_1 = person.segment(1).a(1);
b_1 = person.segment(1).b(1);
h = neck_height;
k = 0.30*c;

if k>h
    k=h;
end
%from fortran code

%% Calculations

% Volume
v_e = 4.66493*a*b*c;
v_c = PI*a_1*b_1*h;
v_p = PI/2*b*c*a_1*(PI/2-b_1*(1-k/c)/b-asin(1-k/c));
volume = v_e + v_c - v_p;

% Mass
m_e = gamma_e*v_e;
m_c = gamma_c*v_c;
m_p = gamma_e*v_p;
mass = m_e + m_c - m_p;

% Mass centroid:
xc = 0;
yc = 0;
zc = (m_e*(c+h-k)+m_c*h/2-m_p*(h-k/3))/mass;

% Moments of inertia:
I_xe = m_e*(0.19473*b^2+0.2351*c^2);    %fortran code has 0.2351, Hatze 79 has 0.23511
I_ye = m_e*(0.221*a^2+0.2351*c^2);  %fortran code has 0.221, Hatze 79 has 0.211
I_ze = m_e*(0.221*a^2+0.19473*b^2);
I_xc = m_c*(3*b_1^2+h^2)/12;
I_yc = m_c*(3*a_1^2+h^2)/12;
I_zc = m_c*(a_1^2+b_1^2)/4;

% principal moments of inertia w.r.t centroid
Ip_x = I_xe+m_e*(c+h-k-zc)^2+I_xc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_y = I_ye+m_e*(c+h-k-zc)^2+I_yc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_z = I_ze+I_zc-m_p*(a_1^2+b_1^2)/6;

% principal moments of inertia w.r.t local systems origin
PIOX=Ip_x+mass*zc^2;
PIOY=Ip_y+mass*zc^2;
PIOZ=Ip_z;

%coordinates of origin of axes (origin of distal segment relative to local
%coordinate system of proximal segment)
OX=0;
OY=l*sin(person.segment(1).theta);
OZ=l*cos(person.segment(1).theta);

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot

  plot_elliptic_plate(P,[a_1 b_1],h,'rotate',R,...
    'colour',person.segment(S).colour,...
    'opacity',person.segment(S).opacity(1),...
    'edgeopacity',person.segment(S).opacity(2));

  zf = @(x,y) c*sqrt(1-(y/b).^2).*(1-(x./a).^8);

  N = 21;
  t = linspace(-PI,PI,N);
  d = linspace(0,1,N);

  [tt,dd] = meshgrid(t,d);

  ho = h + c/2 + k;

  x = dd*a.*cos(tt);
  y = dd*b.*sin(tt);
  z = zf(x,y);

  opt = {'facecolor',person.segment(S).colour,...
    'facealpha',person.segment(S).opacity(1),...
    'edgealpha',person.segment(S).opacity(2)};

  
  surf_rotate(x,y,ho+z,P,R,opt{:})
  surf_rotate(x,y,ho-z,P,R,opt{:})

  x = a*sin(t);
  y = b*cos(t);
  z = zf(x,y);
  surf_rotate([x;x],[y;y],ho+[z;-z],P,R,opt{:})

end

end
