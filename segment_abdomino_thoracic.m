function person = segment_abdomino_thoracic(person,S)
%% abdomino-thoracic

person.segment(S).origin = person.q(1:3);
O1 = person.segment(S).origin + person.segment(S).offset;
i_m = person.sex;

N  = 10; % number of disks
Nt =  7; % number of disks for thoracic region

% Indices for thoracic and abdominal (resp.) groups of disks:
indt = 1:Nt;
inda = (Nt+1):N;

% Measurements and variables:
a  = nan(1,N);
b  = nan(1,N);
w  = nan(1,N);

% coefficients for lungs calcs: (how to generalise?)
c      = nan(1,Nt);
c(1:2) = 0.8333;
c(3:6) = 0.2;
c(7)   = 0.4111;

% Densities:
gamma_t = person.density.thoracic_wall(i_m);
gamma_a = person.density.abdomen(i_m);
gamma_l = person.density.lungs(i_m);
gamma_b = person.density.thoracic_wall(i_m); % (? see A2.94)

%% Measurements
%
% Trying to use Hatze's measurements where possible

l = person.meas{1}.length;

d_11 = person.meas{11}.all(21); % AP distance between centre of hip joint & Symphysion
                    % maybe abdomino-pelvic measurement #21

z_h = mean([person.meas{3}.all(4), person.meas{7}.all(4)]); % height between shoulder and O1

% invented: (can't figure which they'd be)
d = 0.140;    % nipple-to-nipple distance
h = 0.55*l;   % height below C5 of nipple
r = 0.060;    % radius of breast

% thorax ML widths (7) and AP thicknesses (10)
X1 = person.meas{S}.widths;
Y1 = person.meas{S}.depths;
person.meas{S}.length = l;

%% Implicit measurements

RR = person.segment(S).Rglobal;

person.segment(S+1).origin = O1+person.segment(S).Rlocal*[0;0;l];

g = (1+0.3*i_m)*d_11;
jj = round(h/(l/N));

% symmetric chest until the end of the lungs:
w(indt)  = Y1(indt)/2;
b(indt)  = Y1(indt)/2;

% interpolate width minus shoulder; implies 10 disks:
a([1, 5:10]) = X1([1, 5:10])/2;
a(4) = a(5);
ii = [2, 3];
a(ii) = a(5)+(0.42*a(5)-a(1))*l.*(4-ii)/N/(0.35*l-z_h)+(2*a(1)-1.42*a(5))*((1/N*l*(4-ii))/(0.35*l-z_h)).^2;

% interpolate asymmetric belly thicknesses:
w(inda) = interp1([Nt N],[Y1(Nt)/2 g],inda);
b(inda) = Y1(inda) - w(inda);

person.segment(S).a = a;
person.segment(S).b = b;

% Lungs:
a2 = a(indt).*(c(1)-c);
b2 = (b(indt)-a(indt)/6).*sqrt(1-(c/c(1)).^2);


%% Calculations

v_e = pi*a(indt).*b(indt)*l/N; % volume of each thoracic disk
v_p = 8/3*a2.*b2*l/N; % volume of lungs in each disk

m_e = gamma_t*v_e; % mass of thorax as if it were solid
m_p = (gamma_t-gamma_l)*v_p; % mass difference between thorax & lungs
m_t = (v_e-v_p)*gamma_t; % mass of thoracic volume without lungs
m_g = v_p*gamma_l; % mass of lungs only

v_1 = pi*a(inda).*w(inda)*l/2/N;
v_2 = pi*a(inda).*b(inda)*l/2/N;
m_1 = gamma_a*v_1;
m_2 = gamma_a*v_2;

v_f = (1-i_m)*4/3*pi*r^3; % breasts (2 hemispheres)
m_f = gamma_b*v_f;

volume = sum(v_e)+sum(v_1+v_2)+v_f;
mass = sum(m_t+m_g)+sum(m_1+m_2)+m_f;

% Mass centroid:
xc = 0;
yc = ( ...
  sum( (m_1+m_2).*0.424.*(b(inda).^2+w(inda).^2)./Y1(inda) ) ...
  + m_f*(b(jj)+3/8*r) ...
  )/mass;
zc = ( ...
    sum( (m_t+m_g).*(21-2*indt)*l/2/N ) ...
  + sum( (m_1+m_2).*(21-2*inda)*l/2/N ) ...
  + m_f*(l-h) ...
  )/mass;

% Moments of inertia:
s = l^2/1200;
I_x = m_e.*(b(indt).^2+s) ...
      - m_p.*(b2.^2/5+s) ...
      + (m_e-m_p).*(yc^2+(l*(21-2*indt)/20-zc).^2);
Ip_x = sum(I_x) ...
       + sum(...
           m_1.*(0.07*w(inda).^2+s+(-0.424*w(inda)-yc).^2+(l*(21-2*inda)./20-zc).^2) ...
         + m_2.*(0.07*b(inda).^2+s+(+0.424*b(inda)-yc).^2+(l*(21-2*inda)./20-zc).^2) ...
       ) ...
       + m_f*(0.2594*r^2+(l-h-zc)^2+(b(jj)+3*r/8-yc)^2);
I2_x = Ip_x;

I_y = m_e.*(a(indt).^2+s) ...
      - m_p.*(0.06857*a2.^2+s+(c.*a(indt)+0.4*a2).^2) ...
      + (m_e-m_p).*(l*(21-2*indt)/20-zc).^2;

I2_y = sum(I_y) + sum( ...
         (m_1+m_2).*((a(inda).^2)/4+s+(l*(21-2*inda)./20-zc).^2) ...
       ) + m_f*(0.4*r^2+(l-h+zc)^2+(d/2)^2);

I_z = m_e.*(a(indt).^2+b(indt).^2/4) ...
      - m_p.*(0.06857*a2.^2+(b2.^2)/5+(c.*a(indt)+0.4*a2).^2) ...
      + (m_e-m_p)*yc^2;

I2_z = sum(I_z) + sum( ...
         m_1.*(0.07*w(inda).^2+(a(inda).^2)/4+(-0.424*w(inda)-yc).^2) ...
         + m_2.*(0.07*b(inda).^2+(a(inda).^2)/4+(0.424*b(inda)-yc).^2) ...
       ) + m_f*(0.2594*r^2+(b(jj)+3*r/8-yc)^2+(d/2)^2);

I_yz = (m_e-m_p).*(-yc).*(l*(21-2*indt)/20-zc);

I2_yz = sum(I_yz) + sum( ...
    (l*(21-2*inda)./20-zc).*(...
       m_1.*(-0.424*w(inda)-yc) ...
      +m_2.*(+0.424*b(inda)-yc)) ...
  ) + m_f*(b(jj)+3*r/8-yc)*(l-h-zc);

Ip_y = (I2_y+I2_z)/2+sqrt(((I2_y-I2_z)^2)/4+I2_yz^2);
Ip_z = (I2_y+I2_z)/2-sqrt(((I2_y-I2_z)^2)/4+I2_yz^2);

theta = atan(I2_yz/(I2_y-Ip_z));

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).theta = theta;
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot

  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};
  optl = {'opacity',min(1,2*person.segment(S).opacity(1)),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  % thorax:
  for ii = indt
    ph = l-ii*l/N; % plate height

    plot_elliptic_plate(O1+RR*[0;0;ph],[a(ii) w(ii)],l/N,opt{:},'rotate',RR);

    % lungs:
    plot_parabolic_plate(O1+RR*[ c(ii)*a(ii);0;ph],[ a2(ii) b2(ii)],l/N,optl{:});
    plot_parabolic_plate(O1+RR*[-c(ii)*a(ii);0;ph],[-a2(ii) b2(ii)],l/N,optl{:});
  end

  % abdomen:
  for ii = inda
    ph = l-ii*l/N; % plate height
    plot_elliptic_plate(O1+RR*[0;0;ph],[a(ii) -w(ii)],l/N,'segment',[0 0.5],opt{:},'rotate',RR)
    plot_elliptic_plate(O1+RR*[0;0;ph],[a(ii)  b(ii)],l/N,'segment',[0 0.5],opt{:},'rotate',RR)
  end

  % breasts:
  if i_m == 0 % female
    plot_sphere(O1+RR*[+d/2; b(jj); l-h],r,'latrange',[-1 1],'N',[20 10],opt{:})
    plot_sphere(O1+RR*[-d/2; b(jj); l-h],r,'latrange',[-1 1],'N',[20 10],opt{:})
  end

end
