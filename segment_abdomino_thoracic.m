function person = abdomino_thoracic(person,S)
%% abdomino-thoracic

O1 = person.origin{S} + person.offset{S};
i_m = person.sex;

N  = 10; % number of disks
Nt =  7; % number of disks for thoracic region

% Indices for thoracic and abdominal (resp.) groups of disks:
indt = 1:Nt;
inda = (Nt+1):N;

% Measurements and variables:
Y1 = nan(1,N);
X1 = nan(1,N);
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

l = person.meas{1}.all(20);     % length of abdomino-thorasic section; has to be this
                    % value as it's the largest in the set

d_11 = person.meas{11}.all(21); % AP distance between centre of hip joint & Symphysion
                    % maybe abdomino-pelvic measurement #21

z_h = mean([person.meas{3}.all(4), person.meas{7}.all(4)]); % height between shoulder and O1

% invented: (can't figure which they'd be, if any)
d = 0.140;    % nipple-to-nipple distance
h = 0.55*l;   % height below C5 of nipple
r = 0.060;    % radius of breast

% thorax ML widths (7) and AP thicknesses (10)
X1 = person.meas{1}.widths;
Y1 = person.meas{1}.depths;
person.meas{1}.length = l;

%% Implicit measurements

O2 = O1+[0;0;l];
person.origin{2} = O2;

g = (1+0.3*i_m)*d_11;
jj = round(h/(l/N));

% symmetric chest until the end of the lungs:
w(indt)  = Y1(indt)/2;
b(indt)  = Y1(indt)/2;

% interpolate width minus shoulder; implies 10 disks:
a([1 5:10]) = X1([1 5:10])/2;
a(4) = a(5);
ii = [2 3];
a(ii) = a(5)+(0.42*a(5)-a(1))*l.*(4-ii)/N/(0.35*l-z_h)+(2*a(1)-1.42*a(5))*((1/N*l*(4-ii))/(0.35*l-z_h)).^2;

% interpolate asymmetric belly thicknesses:
w(inda) = interp1([Nt N],[Y1(Nt)/2 g],inda);
b(inda) = Y1(inda) - w(inda);

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

V = sum(v_e)+sum(v_1+v_2)+v_f;
M = sum(m_t+m_g)+sum(m_1+m_2)+m_f;

% Mass centroid:
xc = O1(1);
yc = O1(2)+( ...
  sum( (m_1+m_2).*0.424.*(b(inda).^2+w(inda).^2)./Y1(inda) ) ...
  + m_f*(b(jj)+3/8*r) ...
  )/M;
zc = O1(3)+( ...
    sum( (m_t+m_g).*(21-2*indt)*l/2/N ) ...
  + sum( (m_1+m_2).*(21-2*inda)*l/2/N ) ...
  + m_f*(l-h) ...
  )/M;

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

disp('-------------------------')
disp('Abdomino-thoracic section')
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',M)
fprintf('Volume:   %1.4f m^3\n',V)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Theta:    %2.2f°\n',theta*180/pi)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [M,V,1000*xc,1000*yc,1000*zc,theta,Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot

opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
optl = {'opacity',min(1,2*person.opacity{S}(1)),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};

% thorax:
for ii = indt
  ph = l-ii*l/N; % plate height
  
  plot_elliptic_plate(O1+[0;0;ph],[a(ii) w(ii)],l/N,opt{:});
  
  % lungs:
  plot_parabolic_plate(O1+[ c(ii)*a(ii);0;ph],[ a2(ii) b2(ii)],l/N,optl{:});
  plot_parabolic_plate(O1+[-c(ii)*a(ii);0;ph],[-a2(ii) b2(ii)],l/N,optl{:});
end

% abdomen:
for ii = inda
  ph = l-ii*l/N; % plate height
  plot_elliptic_plate(O1+[0;0;ph],[a(ii) -w(ii)],l/N,'segment',[0 0.5],opt{:})
  plot_elliptic_plate(O1+[0;0;ph],[a(ii)  b(ii)],l/N,'segment',[0 0.5],opt{:})
end

% breasts:
if i_m == 0 % female
  plot_sphere(O1+[+d/2; b(jj); l-h],r,'latrange',[-1 1],opt{:})
  plot_sphere(O1+[-d/2; b(jj); l-h],r,'latrange',[-1 1],opt{:})
end

% principle axes:
%plot_coord(O1,'index','1','axes','yz','length',0.1,'rotate',[theta*180/pi 0 0])

% centroid:
%plot3(xc,yc,zc,'.k', 'markersize',20)

end
