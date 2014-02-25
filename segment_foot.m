function person = segment_foot(person,S)

P = person.segment(S).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;
i_m = person.sex;

%% Foot

N = 100; % number of disks
ind = 1:N;

%% Measurements

b  = person.meas{S}.all(1);  %width at toes
c  = person.meas{S}.all(2);     %width at heel
h1 = person.meas{S}.all(3);     %height of big toe
h2 = person.meas{S}.all(4) - h1;    %height of ankle
L  = person.meas{S}.all(5);     %length of foot
d  = person.meas{S}.all(6); % used for ground contact point 

a  = person.meas{S-1}.diam(end); % last leg width measurement

%% Calculations

l_i = a+(L-a)*ind./N;       %A2.106
b_i = a+(b-a)*ind./N;
c_i = a+(c-a)*ind./N;

v    = l_i.*(h2/N).*(b_i+c_i)/2;
v_11 = L*h1*(5*b+c)/36;
v_12 = L*h1*(b+5*c)/36;
v_13 = L*h1*(b+c)/4;
v_14 = 0.5*h2*(L*(a+b+c)/3+a^2);  %from the fortran code - subroutine SEG147, equation S4V14
volume = v_11 + v_12 + v_13 + v_14;

m    = person.density.foot(ind).*v;
m_11 = person.density.heel*v_11;
m_12 = person.density.heel*v_12;
m_13 = person.density.sole*v_13;
mass = m_11 + m_12 + m_13 + sum(m);

% Mass centroid:
xc = 0;
ycm = (m_11+m_12)*(-h2-h1*0.75) +...
    m_13*(-h2-h1*0.25) + ...
    sum(m.*(-ind.*h2/N));
zcm = ...
  m_11*L*(-2/3+(7*b+2*c)/(45*b+9*c)) + ...
  m_12*L*(b+8*c)/(9*b+45*c)         + ...
  m_13*L*(-2/3+(b+2*c)/(3*b+3*c))   + ...
  sum(m.*(-a/2 - ind.*(2*L/3-a/2)/100 + ...
       l_i.*(b_i+2*c_i)./(3*b_i+3*c_i)) ...
  );
yc = ycm/mass;
zc = zcm/mass;

% Moments of inertia:
%%ALL REMAINING 'ERRORS' IN THE MOMENTS OF INERTIA RESULT 
%%FROM HATZE'S USE OF 0.667 IN HIS CODE INSTEAD OF 2/3

Ix_i = m.*(l_i).^2/18.*(b_i.^2+4*b_i.*c_i+c_i.^2)./((b_i+c_i).^2);   
Iz_i = m.*(b_i.^2+c_i.^2)/24;
Iy_i = Ix_i + Iz_i;

% Principal moments of inertia w.r.t. centroid; 
bb=(b+2*c)/3;
cc=(2*b+c)/3;

Ip_x = m_11*( ...
    (L/3)^2*(b^2+4*b*cc+cc^2)/(18*(b+cc)^2) + ...
    (h2+3*h1/4+yc)^2 + (2*L/3 - L*(7*b+2*c)/(45*b+9*c) + zc )^2 ...
  ) + m_12*( ...
    (L/3)^2*(bb^2+4*bb*c+c^2)/(18*(bb+c)^2) + ...
    (h2+3*h1/4+yc)^2 + (L*(b+8*c)/(9*b+45*c) - zc)^2 ...
  ) + m_13*( ...
    (L^2)*(b^2+4*b*c+c^2)/(18*(b+c)^2)+...
    (h2+h1/4+yc)^2 + (2*L/3-L*(b+2*c)/(3*(b+c))+zc)^2 ...
  ) + sum( ...
    Ix_i + m.*((ind.*h2/N+yc).^2 + ...
    (a/2 + (2*L/3-a/2).*ind./N - l_i.*(b_i+2*c_i)./(3.*(b_i+c_i)) + zc).^2) ...  %extra end bracket req in eq A2.108
  );
 
I_y1 = ...
  m_11*(...
    (L/3)^2*(b^2+4*b*cc+(cc)^2)/(18*(b+cc)^2) + ...
    (b^2+cc^2)/24+(2*L/3-L*(7*b+2*c)/(9*(5*b+c))+zc)^2 ...
  ) + m_12*( ...
    (L/3)^2*((bb)^2 + 4*c*(bb)+c^2)/(18*(bb+c)^2) + ...
    (c^2+bb^2)/24+(L*(b+8*c)/(9*(b+5*c))-zc)^2 ...
  ) + m_13*( ...
    L^2*(b^2+4*b*c+c^2)/(18*(b+c)^2) + (b^2+c^2)/24 + (2*L/3-L*(bb/(b+c))+zc)^2 ...
  ) + sum( ...
    Iy_i + m.*((a/2 + (2*L/3-a/2).*ind./N - l_i.*(b_i+2*c_i)./(3.*(b_i+c_i)) + zc).^2) ...
  );

I_z1 = ...
  m_11*( (b^2+cc^2)/24 + (h2+3*h1/4+yc)^2 ) + ...
  m_12*( (bb^2+c^2)/24 + (h2+3*h1/4+yc)^2 ) + ...
  m_13*( (b^2+c^2)/24 + (h2+h1/4+yc)^2) + ...
  sum( Iz_i+m.*(ind.*h2/N+yc).^2 );

I_y1z1 = ...
  m_11*(-h2-3*h1/4-yc)*(-2*L/3+L*(7*b+2*c)/(9*(5*b+c))-zc) + ...
  m_12*(-h2-3*h1/4-yc)*(L*(b+8*c)/(9*b+45*c)-zc) + ...
  m_13*(-h2-h1/4-yc)*(-2*L/3+L*(b+2*c)/(3*(b+c))-zc) + ...
  sum( m.*(-ind.*h2/N-yc).* ...
      (-a/2 - (2*L/3-a/2).*ind./N + l_i.*(b_i+2*c_i)./(3*(b_i+c_i))-zc) ...
  );

Ip_y = (I_y1+I_z1)/2+sqrt(1/4*(I_y1-I_z1)^2+I_y1z1^2);
Ip_z = (I_y1+I_z1)/2-sqrt(1/4*(I_y1-I_z1)^2+I_y1z1^2);
theta = atan(I_y1z1/(I_z1-Ip_y)); 

%centroid w.r.t local coordinate systems (since principal axes differ from
%original segment axes)
xbc = 0;
ybc = yc*cos(theta)+zc*sin(theta);
zbc = zc*cos(theta)-yc*sin(theta);

%principal moments of inertia w.r.t local systems origin
PIOX = Ip_x+mass*(ybc^2+zbc^2);
PIOY = Ip_y+mass*(xbc^2+zbc^2);
PIOZ = Ip_z+mass*(xbc^2+ybc^2);

%coordinates of ground contact points (Heel and Toe)
VXH = 0;
r = person.meas{S-1}.ankle/2;
VYH = -(h2+h1)*cos(theta) + (0.208*L+r/2)*sin(theta);
VZH =  (0.208*L+r/2)*cos(theta) + (h2+h1)*sin(theta);
VXT = 0;
VYT = -(h2+h1)*cos(theta) + (-d+r/2)*sin(theta);
VZT =  (-d+r/2)*cos(theta) + (h2+h1)*sin(theta);


person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];
person.segment(S).theta = theta;

%% Plot

if person.plot || person.segment(S).plot

  opt1  = {'opacity',person.segment(S).opacity(1),'edgeopacity',0,'colour',person.segment(S).colour};
  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  for ii = ind
    ph = -ii*h2/100; % plate height
    s = -a/2-ii/100*(2/3*L-a/2);
    plot_trapzoidal_plate(P+R*[0;ph;s+l_i(ii)/2],c_i(ii)/2,b_i(ii)/2,l_i(ii),h2/100,opt1{:},'rotate',R)
  end

  plot_trapzoidal_plate(P+R*[0;-h2;-L/6],c/2,b/2,L,-h1/2,opt{:},'rotate',R)
  plot_trapzoidal_plate(P+R*[0;-h2-h1;L/6],c/2,(c+1/3*(b-c))/2,L/3,h1/2,opt{:},'rotate',R)
  plot_trapzoidal_plate(P+R*[0;-h2-h1;-L/2],(b+1/3*(c-b))/2,b/2,L/3,h1/2,opt{:},'rotate',R)

end

end

