function [calcs, O13] = thigh(O12,i_m,lr,left_thigh_diameters,left_thigh_perimeters,left_thigh_length,left_thigh_length2)

%% Thigh

N = 10; 
indf = 1:N;

%% Densities

wp= 1;
wt= 1;
v = wp/wt-1; % v is the subcutaneous fat indicator;
gamma_i1 = @(ii,i_m) 1000+(30+10*(ii-2))/((1+2*v)^2)+20*i_m; % for i = 1,2,3;
gamma_i2 = @(i_m) 1030+10*i_m; % for i = 4,5,6,7,8,9;
gamma_i3 = @(i_m) 1490+10*i_m; % for i =10;
gamma_0 = 1020 + 30/((1+2*v)^2)+20*i_m; 

gamma_i = @(i_m) [gamma_i1(1:3,i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i3(i_m)];

%% Measurements

a = left_thigh_diameters/2; % ai 
u = left_thigh_perimeters; 
b = sqrt(((u/pi).^2)/2-a.^2); % b
l = left_thigh_length;
l_1 = left_thigh_length2; 

%% Calculations

% Mass 
v_i = pi*a.*b*l_1/N; % volume of each forearm disk
m_i = gamma_i(i_m).*v_i; % mass of each forearm disk
a_1 = a(1);
b_1 = b(1);
v_0 = 2*pi*a_1*b_1*(l-l_1)/3;
m_0 = gamma_0*v_0;

v=sum(v_i)+v_0;
m=sum(m_i)+m_0;

% Mass centroid:
xc = O12(1);
yc = O12(2);
zc = O12(3) - (m_0*0.4*(l-l_1)+sum(l-l_1+l_1*(2*indf-1)/20))/m;

% Moments of inertia:
I_x0 = m_0*(b_1^2/4+0.0686*l_1^2);
I_y0 = m_0*(0.15*a_1^2+0.0686*l_1^2);
I_z0 = m_0*(0.15*a_1^2+b_1^2/4);
I_xi = m_i.*(3*b.^2+(l/10)^2)/12;
I_yi = m_i.*(3*a.^2+(l/10)^2)/12;
I_zi = m_i.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = I_x0+m_0.*(-0.4*(l-l_1)-zc).^2+sum(I_xi+m_i.*(l-l_1*(2*indf-1)/20-zc).^2);
Ip_y = I_y0+m_0.*(-0.4*(l-l_1)-zc).^2+sum(I_yi+m_i.*(l-l_1*(2*indf-1)/20-zc).^2);
Ip_z= I_z0+sum(I_zi);

disp('-------------------------')
if lr == 'l'
  disp('Left thigh section')
elseif lr == 'r'
  disp('Right thigh section')
end
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',m)
fprintf('Volume:   %1.4f m^3\n',v)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [m,v,xc,yc,zc, Ip_x,Ip_y,Ip_z];

%% Plot
 
O13 = O12+[0;0;-l];
 
opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};
 
for ii = indf
  ph = l_1-ii*l_1/N; % plate height
  plot_elliptic_plate(O13+[0;0;ph],[a(ii) b(ii)],l_1/N,opt{:})
end

end
