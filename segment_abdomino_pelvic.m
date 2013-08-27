function person = segment_abdomino_pelvic(person,S)

O1 = person.origin{1} + person.offset{S};
person.origin{S} = O1;
i_m = person.sex;
nu  = person.nu;

%% Abdominal pelvic region

Nt = 10;
Na = 7;
Np = 3;
ind_ae = 1:Na;
ind_pe = 1:Np;
ind_p = 4:7;
ind_pt = (Np+1):Nt;
ind_as = (Na+1):Nt;

%% Measurements

% WIDTHS  287 301 318 331 341 344 359 382 326 240
%  PERIM  767 797 844 887 940 954 975
%   MEAS  097 081 228 053

a = person.meas{S}.diam([1:8 8 8])/2; % half-widths. what's with 9 & 10 ??
u = person.meas{S}.perim;  % perimeters

a_h  = a(10);  % ? half hip width (trochanter major bones)
B    = person.meas{S}.all(17+4);  % ? bum depth
d_11 = person.meas{S}.all(17+4);  % ? NB also used in ab-thor
e_11 = person.meas{S}.all(17+3);  % ?

h_hoof = Np/Nt*person.meas{11}.length;

l    = person.meas{S}.length;  % "height" of AP section
c = 0.44; % A2.76
h = l/Nt; % height of each plate

%% Densities

gamma_o  = person.density.penis;
gamma_ot = person.density.thigh_head(i_m,nu); 
gamma_1  = person.density.lower_back(i_m);
gamma_2  = person.density.posterior(i_m);
gamma_3  = person.density.stomach(i_m);
gamma_4  = person.density.pelvis(i_m);
gamma_5  = person.density.buttocks(i_m,nu);

%% Thighs

atl = person.meas{12}.diam/2; 
utl = person.meas{12}.perim; 
btl = sqrt(((utl/pi).^2)/2-atl.^2);
atr = person.meas{15}.diam/2; 
utr = person.meas{15}.perim; 
btr = sqrt(((utr/pi).^2)/2-atr.^2);

O12 = O1 + [-btl(1); 0; -l + h_hoof];
O15 = O1 + [+btr(1); 0; -l + h_hoof];

person.origin{12} = O12;
person.origin{15} = O15;

%% Geometry

g = (1+0.3*i_m)*d_11;

f_1(ind_pt) = c*a_h*(1+sqrt(1-(0.88558-0.22883*(ind_pt-4)).^2));
r = sqrt(e_11^2-(0.037*l)^2)-B-g;

% anterior depth of stomach above the buttocks:
b(ind_pe) = sqrt( 2*(u(ind_pe)/pi-sqrt(0.5*(a(ind_pe).^2+g^2))).^2 - a(ind_pe).^2 );

vv = 1-(0.88558-0.22883*(ind_p-4)).^2;
aa(ind_p) = c*a_h*sqrt(vv);
bb(ind_p) = B*vv;

s_l(ind_p) = sqrt(g^2+(a(ind_p)-f_1(ind_p)).^2);
rr = 2*bb(ind_p)./aa(ind_p);
s_p(ind_p) = 0.5*aa(ind_p).*(...
   sqrt(1+rr.^2) + 1./rr.*log(rr+sqrt(1+rr.^2))...
  );

b(ind_p) = sqrt(0.20264*(...
    u(ind_p) - 2*( c*a_h + s_p(ind_p) + s_l(ind_p) )...
  ).^2 - a(ind_p).^2 );

% correction if necessary for really bad inputs
for ii = 1:length(b)
  if ~isreal(b(ii))
    [ii b(ii)]
    b(ii) = b(ii-1);
  end
end


%% Mass

% m_p = gamma_5*pi/2*B*c*a_h*0.473*l; % two elliptic paraboloids;
% m_ei1 = gamma_1*pi/2*g*a_i*l/10; % i=1,2,3; three semi-elliptical plates;
% m_ti = gamma_2*g*l/10*(a_i+f_1); % i=4,5,6,7,8,9,10; seven trapezoidal plates on the posterior side of the segment;
% m_ei = gamma_3*pi/2*a_i*b_i*l/10; % i=1,2,3,4,5,6,7; seven semi-elliptical plates;
% m_l = gamma_4*0.3*l*(2*a_8*r-(2-pi/2)*a_1t*b_1t); % three special shape plates on the anterior side;
% m_otl = gamma_ot*2*pi*a*b*h/3; % the removed superior parts of left thigh;
% m_otr = gamma_ot*2*pi*a*b*h/3; % the removed superior parts of right thigh;
% 
% m_o = 0.007*i_m; % when A is no lager than 12;
% m_o = gamma_o*i_m*2*pi*(0.005*A-0.045)^3/3; % when A is between 12 and 19;
% m_o = 0.26*i_m; % when A is lager than 19;
% 
% nu = m_p*2/gamma_5 + sum(m_ei1)/gamma_1 + sum(m_ti)/gamma_2 + sum(m_ei)/gamma_3 + m_l*3/gamma_4 + m_o - m_otl - m_otr;
% m = m_p*2 + sum(m_ei1) + sum(m_ti) + sum(m_ei) + m_l*3 + m_o - m_otl - m_otr;
% 
% % Mass centroid:
% xc = O11(0);
% y_p = -B/3-g;
% y_1 = -4*g/3/pi;
% y_2i = -(g/3)*(2*f_l+a_i)/(i_l+a_i);
% y_3i = 4*b_i/3/pi;
% y_4 = (2*(a_i(8)-a_lt)*b_lt*(r-b_lt/2)+a_i(8)*(r-b_lt)^2+1.571*a_lt*b_lt*(r-0.576*b_lt))/(2*a_i(8)*r-(2-pi/2)*a_lt*b_lt);
% y_ot = r-b_lt;
% z_p = -0.737*l;
% z_li = -l*(indf/10-0.05);
% z_4 = -0.85*l;
% z_ot = -0.7*l-0.6*(l_t-l_lt);
% yc = (2*m_p*y_p+y_1*sum(m_ei)+y_2i*sum(m_ti) + y_3i*sum(m_ei) + y_4*m_l*3 + m_o*(r+0.02) - m_otl*z_ot - m_otr*z_ot)/m;
% zc = (m_p*2*z_p + sum(m_ei1)*z_li + sum(m_ti)*z_li + sum(m_ei)*z_li + m_l*3*z_li + m_o - m_otl*z_ot - m_otr*z_ot)/m;
% 
% 
% disp('-------------------------')
% disp('Abdomino pelvic section')
% disp('-------------------------')
% fprintf('Mass:     %2.3f kg\n',sum(m))
% fprintf('Volume:   %1.4f m^3\n',sum(nu))
% fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
% fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

person.segment(S).theta = 10; % needs to be calculated

%% Plot

if person.plot
  
  opt = {'colour',person.color{S},'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2)};
  
  % buttocks left and right
  
  plot_elliptic_paraboloid(O1+[-c*a_h;-g;-l+h_hoof-0.037*l],0.437*l,B,'rotate',[90 0 0],'N',[20 7],opt{:})
  plot_elliptic_paraboloid(O1+[+c*a_h;-g;-l+h_hoof-0.037*l],0.437*l,B,'rotate',[90 0 0],'N',[20 7],opt{:})
  
  % posterior: 3 semi-elliptical plates
  for ii = ind_pe
    plot_elliptic_plate(O1+[0;0;-ii*h],[a(ii) g],h,'segment',[0.5 1],opt{:})
  end
  
  % posterior: 7 trapezoidal plates
  for ii = ind_pt
    plot_trapzoidal_plate(O1+[0;-g/2;-ii*h],2*a(ii),2*f_1(ii),g,h,opt{:})
  end
  
  % anterior: 7 semi-elliptical plates
  for ii = ind_ae
    plot_elliptic_plate(O1+[0;0;-ii*h],[a(ii) b(ii)],h,'segment',[0 0.5],opt{:})
  end
  
  % anterior: 3 "special shape" plates
  for ii = ind_as
    plot_special_plate(O1+[0;0;-ii*h],2*a(8),r,atl(1),btl(1),h,opt{:})
  end
  
end
