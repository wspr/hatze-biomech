%% abdomino-thoracic

clear all

N  = 10; % number of disks
Nt =  7; % number of disks for thoracic region

male   = 1;
female = 0;

% Measurements and variables:
Y1 = nan(1,N);
X1 = nan(1,N);
w  = nan(1,N);
a  = nan(1,N);
b  = nan(1,N);
c  = nan(1,Nt);
a2 = nan(1,Nt);
b2 = nan(1,Nt);

% volumes of elements:
v_1 = nan(1,N);
v_2 = nan(1,N);
v_e = nan(1,Nt);

% masses of elements:
m_1 = nan(1,N);
m_2 = nan(1,N);
m_e = nan(1,Nt);
m_p = nan(1,Nt);

% moments of inertia
I_x = nan(1,Nt);
I_y = nan(1,Nt);
I_z = nan(1,Nt);
I_yz = nan(1,Nt);

% coefficients for lungs calcs:
c(1:2) = 0.8333;
c(3:6) = 0.2;
c(7)   = 0.4111;

% Densities:
gamma_t = @(i_m) 1080+60*i_m; % thoracic wall
gamma_a = @(i_m) 1000+30*i_m; % abdomen
gamma_l = 300; % lungs
gamma_b = 1200; % breasts (? see A2.94)

% Indices for thoracic and abdominal (resp.) groups of disks:
indt = 1:Nt;
inda = (Nt+1):N;

%% Measurements
%
% Trying to use Hatze's measurements where possible

i_m = male;

hatze1 = [...
  101 273 275 264 251 245 262 ... 1:7 ML
  178 ... 8   ???
  084 095 117 149 168 192 190 184 174 174 ... 9:18 AP
  170 ... 19  ???
  452 ... 20  l
  213 ... 21  ???
  ]/1000;
hatze2 = [139 184 214 055]/1000;
hatze11(21) = 0.053;

l = hatze1(20);     % length of abdomino-thorasic section; has to be this
                    % value as it's the largest in the set

d_11 = hatze11(21); % AP distance between centre of hip joint & Symphysion
                    % maybe abdomino-pelvic measurement #21

z_h = hatze2(4);    % height between shoulder and O1
                    % guessing it's the smallest of the shoulder
                    % measurements

% invented: (can't figure which they'd be, if any)
d = 0.15;    % nipple-to-nipple distance
h = l/2;     % height below C5 of nipple
r = 0.05;    % radius of breast

ML = 0;
AP = 9;

% thorax ML widths (7)
X1(1)  = hatze1(ML+1);
X1(5)  = hatze1(ML+2);
X1(6)  = hatze1(ML+3);
X1(7)  = hatze1(ML+4);
X1(8)  = hatze1(ML+5);
X1(9)  = hatze1(ML+6);
X1(10) = hatze1(ML+7);

% thorax AP thicknesses (10)
Y1(1)  = hatze1(AP);
Y1(2)  = hatze1(AP+1); 
Y1(3)  = hatze1(AP+2); 
Y1(4)  = hatze1(AP+3);
Y1(5)  = hatze1(AP+4); 
Y1(6)  = hatze1(AP+5); 
Y1(7)  = hatze1(AP+6); 
Y1(8)  = hatze1(AP+7); 
Y1(9)  = hatze1(AP+8); 
Y1(10) = hatze1(AP+9);

%% Implicit measurements

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
a2(indt) = a(indt).*(c(1)-c(indt));
b2(indt) = (b(indt)-a(indt)/6).*sqrt(1-(c(indt)/c(1)).^2);


%% Calculations

v_e(indt) = pi*a(indt).*b(indt)*l/N;
m_e(indt) = gamma_t(i_m)*v_e(indt);
m_p(indt) = (gamma_t(i_m)-gamma_l)*8/3*a2(indt).*b2(indt)*l/N;

v_1(inda) = pi*a(inda).*w(inda)*l/2/N;
v_2(inda) = pi*a(inda).*b(inda)*l/2/N;
m_1(inda) = gamma_a(i_m)*v_1(inda);
m_2(inda) = gamma_a(i_m)*v_2(inda);

v_f = (1-i_m)*4/3*pi*r^3; % breasts (2 hemispheres)
m_f = gamma_b*v_f;

V = sum(v_e(indt))+sum(v_1(inda)+v_2(inda))+v_f;
M = sum(m_e(indt)+m_p(indt))+sum(m_1(inda)+m_2(inda))+m_f;

fprintf('Mass:   %2.3f kg\n',M)
fprintf('Volume: %1.4f m^3\n',V)

% Mass centroid:
xc = 0;
yc = ( ...
  sum( (m_1(inda)+m_2(inda))*0.424.*(b(inda).^2+w(inda).^2)./Y1(inda) ) ...
  + m_f*(b(jj)+3/8*r) ...
  )/M;
zc = ( ...
    sum( (m_e(indt)-m_p(indt)).*(21-2*indt)*l/20 ) ...
  + sum( (m_1(inda)+m_2(inda)).*(21-2*inda)*l/20 ) ...
  + m_f*(l-h) ...
  )/M;

s = l^2/1200;

%% Plot

hfig = figure(1); clf; hold on
set(hfig,'color','white')

opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};

% thorax:
for ii = indt
  ph = l-ii*l/N; % plate height
  
  plot_elliptic_plate([0,0,ph],[a(ii) w(ii)],l/N,opt{:});
  
  % lungs:
  plot_parabolic_plate([ c(ii)*a(ii),0,ph],[ a2(ii) b2(ii)],l/N,optl{:});
  plot_parabolic_plate([-c(ii)*a(ii),0,ph],[-a2(ii) b2(ii)],l/N,optl{:});
end

% abdomen:
for ii = inda
  ph = l-ii*l/N; % plate height
  plot_elliptic_plate([0,0,ph],[a(ii) -w(ii)],l/N,'segment',[0 0.5],opt{:})
  plot_elliptic_plate([0,0,ph],[a(ii)  b(ii)],l/N,'segment',[0 0.5],opt{:})
end

% breasts:
if i_m == female
  plot_sphere([+d/2, b(jj), l-h],r,opt{:})
  plot_sphere([-d/2, b(jj), l-h],r,opt{:})
end

plot3(xc,yc,zc,'.k', 'markersize',40)

axis equal
view(170,15)
axis off
zoom(1.7)