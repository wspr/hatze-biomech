function person = segment_abdomino_pelvic(person,S)

O1 = person.segment(1).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;
person.segment(S).origin = O1;
i_m = person.sex;
nu  = person.nu;
age = person.age;

%% Abdominal pelvic region

Nt = 10;
Na = 7;
Np = 3;
ind = 1:Nt;
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
e_11 = person.meas{S}.all(17+3);  % ? certainly not

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

a_1t = 0.5*( atl(1) + atr(1) );
b_1t = 0.5*( btl(1) + btr(1) );


%% THESE NEED TO BE FIXED: offset should be 0.5*(pelvis_width-thigh_width)
O12 = O1 + person.segment(S).Rglobal*[-btl(1); 0; -l + h_hoof];
O15 = O1 + person.segment(S).Rglobal*[+btr(1); 0; -l + h_hoof];

person.segment(12).origin = O12;
person.segment(15).origin = O15;

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

v_p = pi/2*B*c*a_h*0.473*l; % two elliptic paraboloids;
v_ee(ind_pe) = pi/2*g*a(ind_pe)*l/10; % i=1,2,3; three semi-elliptical plates;
v_t(ind_pt) = g*l/10*(a(ind_pt)+f_1(ind_pt)); % i=4,5,6,7,8,9,10; seven trapezoidal plates on the posterior side of the segment;
v_e(ind_ae) = pi/2*a(ind_ae).*b*l/10; % i=1,2,3,4,5,6,7; seven semi-elliptical plates;
v_l = 0.3*l*(2*a(8)*r-(2-pi/2)*a_1t*b_1t); % three special shape plates on the anterior side;
v_otl = 2*pi*atl(1)*btl(1)*h_hoof/3; % the removed superior parts of left thigh;
v_otr = 2*pi*atr(1)*btr(1)*h_hoof/3; % the removed superior parts of right thigh;

m_p = gamma_5*v_p;
m_ee = gamma_1*v_ee;
m_t = gamma_2*v_t;
m_e = gamma_3*v_e;
m_l = gamma_4*v_l;
m_otl = gamma_ot*v_otl;
m_otr = gamma_ot*v_otr;

% NB: ambiguity between h_hoof for pelvis and the separate h terms for each
% thigh

if age <= 12
  m_o = 0.007*i_m;
elseif age <=19
  m_o = gamma_o*i_m*2*pi*(0.005*age-0.045)^3/3;
else
  m_o = 0.26*i_m;
end

% NB: volume v_o is ambiguous

v = v_p*2 + sum(v_ee) + sum(v_t) + sum(v_e) + v_l*3 + m_o/gamma_o - v_otl - v_otr;
m = m_p*2 + sum(m_ee) + sum(m_t) + sum(m_e) + m_l*3 + m_o - m_otl - m_otr;


% Mass centroid:
xc = 0;
y_p = -B/3-g;
y_1 = -4*g/3/pi;
y_2i(ind_pt) = -(g/3)*(2*f_1(ind_pt)+a(ind_pt))/(f_1(ind_pt)+a(ind_pt));
y_3i = 4*b/3/pi;
y_4 = (2*(a(8)-a_1t)*b_1t*(r-b_1t/2)+a(8)*(r-b_1t)^2+1.571*a_1t*b_1t*(r-0.576*b_1t))/(2*a(8)*r-(2-pi/2)*a_1t*b_1t);
y_ot = r-b_1t;
z_p = -0.737*l;
z_1 = -l*(ind/Nt-0.05); % what is "-0.05" ?
z_4 = -0.85*l;
z_ot = -0.7*l-0.6*h;

ycm = 2*m_p*y_p+y_1*sum(m_e)+sum(y_2i.*m_t) + sum(y_3i.*m_e) + y_4*m_l*3  - (m_otl*y_ot - m_otr*y_ot)+ m_o*(r+0.02);
zcm = 2*m_p*z_p + sum( z_1(ind_pe).*(m_e(ind_pe)+m_ee(ind_pe)) ) + sum(z_1(ind_ae).*m_e(ind_ae)) + sum(z_1(ind_pt).*m_t(ind_pt)) + 3*m_l*z_4 + m_o - (m_otl + m_otr)*z_ot;

% NB: m_e*z_1 seems to be counted TWICE for i=1:3 ?

yc = ycm/m;
zc = zcm/m;

% Moments of inertia:
f_Ti(ind_pt) = (g^2/18)*(a(ind_pt).^2+4*a(ind_pt).*f_1(ind_pt)+f_1(ind_pt).^2)/((a(ind_pt)+f_1(ind_pt)).^2);
g_Ti(ind_pt) = (a(ind_pt).^2+f_1(ind_pt).^2)/6;

Ip_x = ...
  2*m_p*((3*(0.437*l)^2+B^2)/18+(y_p-yc).^2+(z_p-zc).^2) + ...
  sum(m_ee(1:3).*(0.07*g^2+l^2/1200+(y_1-yc).^2+(z_1(1:3)-zc).^2)) + ...
  sum(m_e(1:7).*(0.07*b(1:7).^2+l^2/1200+(y_3i(1:7)-yc).^2+(z_1(1:7)-zc).^2)) + ...
  sum(m_t(4:10).*(f_Ti(4:10)+(y_2i(4:10)-yc).^2+(z_1(4:10)-zc).^2)) + ...
  3*m_l*((r^2+(0.3*l)^2)/12+(y_4-yc).^2+(z_4-zc).^2) + ...
  -(m_otl + m_otr)*(0.25*b_1t^2+0.0686*h^2+(y_ot-yc)^2+(z_ot-zc)^2) + ... % F2.9 the height of the special shape; h=l_t - l_1t;
  m_o*((r+0.02-yc)^2+(-0.7*l-0.05-zc)^2);

I_y = ...
  2*m_p*( ((0.437*l)^2+(c*a_h)^2)/6+(c*a_h)^2+(z_p-zc).^2 ) + ...
  sum(m_ee(1:3).*(0.25*a(1:3).^2+l^2/1200+(z_1(1:3)-zc).^2)) + ...
  sum(m_e(1:7).*(0.25*a(1:7).^2+l^2/1200+(z_1(1:7)-zc).^2)) + ...
  sum(m_t(4:10).*(g_Ti(4:10)+(z_1(4:10)-zc).^2)) + ...
  3*m_l*(((0.3*l)^2+4*a(8)^2)/12+(z_4-zc)^2) + ...
  -(m_otl+m_otr)*(0.15*a_1t^2+0.0686*h ^2+(a(8)-a_1t)^2+(z_ot-zc)^2) + ...
  m_o*(-0.7*l-0.05-zc)^2;

I_z = ...
  2*m_p*((3*(c*a_h)^2+B^2)/18+(y_p-yc).^2+(c*a_h)^2)...
  +sum(m_ee(1:3).*(0.07*g^2+0.25*a(1:3).^2+(y_1-yc).^2))...
  +sum(m_e(1:7).*(0.07*b(1:7).^2+0.25*a(1:7).^2+(y_3i(1:7)-yc).^2))...
  +sum(m_t(4:10).*(f_Ti(4:10)+g_Ti(4:10)+(y_2i(4:10)-yc).^2))...
  +3*m_l*((r^2+4*a(8)^2)/12+(y_4-yc)^2)...
  -(m_otl + m_otr)*(0.15*a_1t^2+0.25*b_1t^2+(a(8)-a_1t)^2+(y_ot-yc)^2)...
  +m_o*(r+0.02-yc)^2;

I_yz = ...
  2*m_p*(y_p-yc)*(z_p-zc)...
  +sum(m_ee(1:3).*(y_1-yc).*(z_1(1:3)-zc))...
  +sum(m_e(1:7).*(y_3i(1:7)-yc).*(z_1(1:7)-zc))...
  +sum(m_t(4:10).*(y_2i(4:10)-yc).*(z_1(4:10)-zc))...
  +3*m_l*(y_4-yc)*(z_4-zc)...
  -(m_otl + m_otr)*(y_ot-yc)*(z_ot-zc)...
  +m_o*(r+0.02-yc)*(-0.7*l-0.05-zc);

Ip_y =(I_y+I_z)/2+sqrt(1/4*(I_y-I_z)^2+I_yz^2);
Ip_z =(I_y+I_z)/2-sqrt(1/4*(I_y-I_z)^2+I_yz^2);

theta = atan(I_yz/(I_z-Ip_y));


person.segment(S).volume = v;
person.segment(S).mass = m;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x, Ip_y, Ip_z];
person.segment(S).theta = theta;


%% Plot

if person.plot || person.segment(S).plot

  opt = {'colour',person.segment(S).colour,'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2)};

  % buttocks left and right

  plot_elliptic_paraboloid(O1+R*[-c*a_h;-g;-l+h_hoof-0.037*l],0.437*l,B,'rotate',R*rotation_matrix_zyx([90 0 0]),'N',[20 7],opt{:})
  plot_elliptic_paraboloid(O1+R*[+c*a_h;-g;-l+h_hoof-0.037*l],0.437*l,B,'rotate',R*rotation_matrix_zyx([90 0 0]),'N',[20 7],opt{:})

  % posterior: 3 semi-elliptical plates
  for ii = ind_pe
    plot_elliptic_plate(O1+R*[0;0;-ii*h],[a(ii) g],h,'segment',[0.5 1],'rotate',R,opt{:})
  end

  % posterior: 7 trapezoidal plates
  for ii = ind_pt
    plot_trapzoidal_plate_xy(O1+R*[0;-g/2;-ii*h],2*a(ii),2*f_1(ii),g,h,'rotate',R,opt{:})
  end

  % anterior: 7 semi-elliptical plates
  for ii = ind_ae
    plot_elliptic_plate(O1+R*[0;0;-ii*h],[a(ii) b(ii)],h,'segment',[0 0.5],opt{:},'rotate',R)
  end

  % anterior: 3 "special shape" plates
  for ii = ind_as
    plot_special_plate(O1+R*[0;0;-ii*h],2*a(8),r,atl(1),btl(1),h,'rotate',R,opt{:})
  end

end
