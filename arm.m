function [calcs, O5] = arm(O4,i_m,lr,arm_diameters,arm_perimeters,arm_length)

%%  Arm

N = 10; 
indu = 1:N;
 
%% Densities

% for i = 1
gamma_i1 = 1060+40*i_m; 
% for i= 2,3,4,5,6,7,8,9:
gamma_i2 = 1058+20*i_m; 
% for i =10;
gamma_i3 = 1080+20*i_m; 
% for i =0;
gamma_0 = gamma_i3; 

gamma_i = [gamma_i1 gamma_i2 gamma_i2 gamma_i2 gamma_i2 gamma_i2 gamma_i2 gamma_i2 gamma_i2 gamma_i3];

%% Measurements

a = arm_diameters/2;
u = arm_perimeters;
b = sqrt(((u/pi).^2)/2-a.^2);
l = arm_length;

%% Calculations

% volume of each forearm disk;
v_i = pi*a.*b*l/N; 
v_0 = 2*pi*((b(1)/2)^3)/3;
v= sum(v_i)+v_0;

% mass of each forearm disk;
m_i = gamma_i.*v_i; 

% mass of the hemisphere;
m_0 = gamma_0*v_0; 

% total mass of left arm;
m = m_0 + sum(m_i); 
 
% Mass centroid:
xc = 0;
yc = 0;
zc = (m_0*0.375*(b(1)/2) - sum(m_i.*(2*indu-1).*l)/2/N)/m;
 
% Moments of inertia:
I_x = m.*(3*(b.^2)+(l/10).^2)/12; 
I_y = m.*(3*(a.^2)+(l/10).^2)/12;
I_z = m.*(3*(b.^2)+(b/10).^2)/12;
 
% principal moments of inertia; 
Ip_x = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_x) +sum(m_i.*(l*(2*indu-1)/20+zc).^2);
Ip_y = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_y) +sum(m_i.*(l*(2*indu-1)/20+zc).^2);
Ip_z= m_0*((b(1)./2)^2)/5 + sum(I_z);
 
disp('-------------------------')
if lr == 'l'
  disp('Left arm section')
elseif lr == 'r'
  disp('Right arm section')
end
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',m)
fprintf('Volume:   %1.4f L\n',v*1000)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [m,v,xc,yc,zc, Ip_x,Ip_y,Ip_z];

%% Plot

O5 = O4+[0;0;-l];

% axes:
plot_coord(O4,'index','4''');
plot_coord(O5,'index','5''');
 
opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};
 
for ii = indu
  ph = l-ii*l/N; % plate height
  plot_elliptic_plate(O5+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
end

% the hemishpere
plot_sphere(O4, b(1)/2, 'longrange',[0 1],opt{:})

end