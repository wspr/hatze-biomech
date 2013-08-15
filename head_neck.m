function [calcs] = head_neck(O2,i_m,hatze1,head_width,head_depth,head_height,neck_height)

%% Head_neck

%% Densities:

gamma_e = 1120;
gamma_c = 1040;

%% Measurements

a = head_width/2;
b = 1.08*head_depth/2;
c = 1.04*head_height/2;
a_1 = hatze1(1)/2;
b_1 = hatze1(9)/2;
h = neck_height;
k = 0.30*c;

%% Calculations

% Mass 
m_e = gamma_e*4.66493*a*b*c;
m_c = gamma_c*pi*a_1*b_1*h;
m_p = gamma_e*4.5708*b*c*a_1*(1.5708-b_1*(c-k)/(c*b)-asin((c-k)/c));
m = m_e + m_c + m_p;

% Mass centroid:
xc = O2(1);
yc = O2(2);
zc = (m_e*(c+h-k)+m_c*h/2-m_p*(h-k/3))/m;

% Moments of inertia:
I_xe = m*(0.19473*b^2+0.23511*c^2);
I_ye = m*(0.211*a^2+0.23511*c^2);
I_ze = m*(0.211*b^2+0.19473*b^2);
I_xc = m*(3*b^2+h^2)/12;
I_yc = m*(3*a^2+h^2)/12;
I_zc = m*(a^2+b^2)/4;

% principal moments of inertia; 
Ip_x = I_xe+m_e*(c+h-k-zc)^2+I_xc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_y = I_ye+m_e*(c+h-k-zc)^2+I_yc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_z = I_ze+I_zc-m_p*(a_1^2+b_1^2)/6;

disp('-------------------------')
disp('head neck section')
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',m)
%fprintf('Volume:   %1.4f m^3\n',v)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [m,xc,yc,zc, Ip_x,Ip_y,Ip_z];

%% Plot

plot_elliptic_plate(O2,[a_1 b_1],h,'colour',[1 0 0],'opacity',0.3,'edgeopacity',0.1);

zf = @(x,y) c*sqrt(1-(y/b).^2).*(1-(x./a).^8);

N = 21;
t = linspace(-pi,pi,N);
d = linspace(0,1,N);

[tt,dd] = meshgrid(t,d);

ho = O2(3) + h + c/2 + k;

x = O2(1)+dd*a.*cos(tt);
y = O2(2)+dd*b.*sin(tt);
z = zf(x,y);

opt = {'facecolor',[1 0 0],'facealpha',0.3,'edgealpha',0.1};

surf(x,y,ho+z,opt{:})
surf(x,y,ho-z,opt{:})

x = O2(1)+a*sin(t);
y = O2(2)+b*cos(t);
z = zf(x,y);
surf([x;x],[y;y],ho+[z;-z],opt{:})

plot3(xc,yc,zc,'k.','markersize',20)