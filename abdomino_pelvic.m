function abdomino_pelvic(O11)

%% Abdomino pelvic
 
clear all

O11 = [0;0;0]; % origin

% For logical selection of gender:
male   = 1;
female = 0;
i_m = male;
N = 10; % number of disks
indf = 1:N; 

% Densities:
gamma_o = 1000;
v = wp/wt-1; % v is the subcutaneous fat indicator;
gamma_ot = 1020 + 30/((1+2*v)^2)+20*i_m; 
gamma_1 = 1090 + 30*i_m;
gamma_2 = 1020+30*i_m;
gamma_3 = 1000 + 40*i_m;
gamma_4 = 1020 + 30*i_m;
gamma_5 = 960+30/(1+4*v)^3+60*i_m;

h =l_t-l_1t;

%% Measurements

hatze11 = [...
  287 301 318 331 341 344 359 382 326 240  097 081 228 053];

perimeters = [767 797 844 887 940 954 975];

%a_h
%l
%u_i
%B
%a_i
%a_1t
%b_1t
%d_11

%% Calculations
c = 0.44;
g = (1+0.3*i_m)*d_11;

f_1 = c*a_h*(1+(1-(0.88558-0-22883*(indf-4))^2)^0.5);
r = (e11^2-(0.037*l)^2)^0.5-B-g;
b_i = (2*(u_i/pi-(0.5*(a_i^2+g^2))^0.5)^2-a_i^2)^0.5; % i=1,2,3;
b_i = (0.20264*(u_i-2*(c*a_h+s_p+s_l))^2-a_i^2)^0.5; % i=4,5,6,7;
s_il = (g^2+(a_i-f_1)^2)^0.5;
s_ip = 0.5*a_i1((1+(2*b_i1/a_i1)^2)^0.5+(a_i1/2*b_i1)*ln(2*b_i1/a_i1+(1+(2*b_i1/a_i1)^2)^0.5));
% a_8 =
a_i1 = c*a_h*(1-(0.88558-0.22883*(indf-4))^2)^0.5;
b_i1 = B*(1-(0.88558-0-22883*(indf-4))^2);

% Mass
m_p = gamma_5*pi/2*B*c*a_h*0.473*l; % two elliptic paraboloids;
m_ei1 = gamma_1*pi/2*g*a_i*l/10; % i=1,2,3; three semi-elliptical plates;
m_ti = gamma_2*g*l/10*(a_i+f_1); % i=4,5,6,7,8,9,10; seven trapezoidal plates on the posterior side of the segment;
m_ei = gamma_3*pi/2*a_i*b_i*l/10; % i=1,2,3,4,5,6,7; seven semi-elliptical plates;
m_l = gamma_4*0.3*l*(2*a_8*r-(2-pi/2)*a_1t*b_1t); % three special shape plates on the anterior side;
m_otl = gamma_ot*2*pi*a*b*h/3; % the removed superior parts of left thigh;
m_otr = gamma_ot*2*pi*a*b*h/3; % the removed superior parts of right thigh;

m_o = 0.007*i_m; % when A is no lager than 12;
m_o = gamma_o*i_m*2*pi*(0.005*A-0.045)^3/3; % when A is between 12 and 19;
m_o = 0.26*i_m; % when A is lager than 19;

v = m_p*2/gamma_5 + sum(m_ei1)/gamma_1 + sum(m_ti)/gamma_2 + sum(m_ei)/gamma_3 + m_l*3/gamma_4 + m_o - m_otl - m_otr;
m = m_p*2 + sum(m_ei1) + sum(m_ti) + sum(m_ei) + m_l*3 + m_o - m_otl - m_otr;

% Mass centroid:
xc = O11(0);
y_p = -B/3-g;
y_1 = -4*g/3/pi;
y_2i = -(g/3)*(2*f_l+a_i)/(i_l+a_i);
y_3i = 4*b_i/3/pi;
y_4 = (2*(a_i(8)-a_lt)*b_lt*(r-b_lt/2)+a_i(8)*(r-b_lt)^2+1.571*a_lt*b_lt*(r-0.576*b_lt))/(2*a_i(8)*r-(2-pi/2)*a_lt*b_lt);
y_ot = r-b_lt;
z_p = -0.737*l;
z_li = -l*(indf/10-0.05);
z_4 = -0.85*l;
z_ot = -0.7*l-0.6*(l_t-l_lt);
yc = (2*m_p*y_p+y_1*sum(m_ei)+y_2i*sum(m_ti) + y_3i*sum(m_ei) + y_4*m_l*3 + m_o*(r+0.02) - m_otl*z_ot - m_otr*z_ot)/m;
zc = (m_p*2*z_p + sum(m_ei1)*z_li + sum(m_ti)*z_li + sum(m_ei)*z_li + m_l*3*z_li + m_o - m_otl*z_ot - m_otr*z_ot)/m;


disp('-------------------------')
disp('Abdomino pelvic section')
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',sum(m))
fprintf('Volume:   %1.4f m^3\n',sum(v))
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

%% Plot

% buttock 1



% buttock 2

% 3 semi-elliptical plates

% 7 trapezoidal plates

% 7 semi-elliptical plates

% 3 "special shape" plates

