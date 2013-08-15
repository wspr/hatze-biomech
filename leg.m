function [calcs, O17] = leg(O16,i_m,lr,left_leg_diameters,left_leg_perimeters,left_leg_length,ankle_size)

%% Leg

N = 10; 
indf = 1:N;

%% Densities

wp= 1;
wt= 1;
v = wp/wt-1; % v is the subcutaneous fat indicator;
gamma_i1 = @(ii,i_m) 1000+(30+10*(ii-2))/((1+2*v)^2)+20*i_m; % for i = 1,2,3;
gamma_i2 = @(i_m) 1034+10*i_m; % for i = 4,5,6,7,8,9;
gamma_i3 = @(i_m) 1490+10*i_m; % for i =10;
gamma_b = 1200;

gamma_i= @(i_m) [gamma_i1(1:3,i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i2(i_m) gamma_i3(i_m)];

%% Measurements

a = left_leg_diameters/2; 
u = left_leg_perimeters;
b = sqrt (((u(indf)/pi).^2)/2-a(indf).^2); 
l = left_leg_length;
r = ankle_size;

%% Calculations
% Mass 
v = pi*a(indf).*b(indf)*l/N; % volume of each forearm disk
m_b = 2*gamma_b*0.92*r^3;
m = gamma_i(i_m).*v+2*m_b; % mass of each forearm disk
mass = sum(m)+m_b;

% Mass centroid:
xc = O16(1);
yc = O16(2);
zc = O16(3) - (-2*m_b*l-sum(m.*l.*(2*indf-1)/20)./mass);

% Moments of inertia:
I_xi = m.*(3*b.^2+(l/10)^2)/12;
I_yi = m.*(3*a.^2+(l/10)^2)/12;
I_zi = m.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = 2*m_b*(0.33*r^2+(l+zc).^2)+sum(I_xi+m.*(l*(2*indf-1)/20+zc).^2);
Ip_y = 2*m_b*(0.1859*r^2+(l+zc).^2+(a(10)+0.196*r^2))+ sum(I_xi+m.*(l*(2*indf-1)/20+zc).^2);
Ip_z = 2*m_b*(0.1859*r^2+(a(10)+0.196*r)^2)+ sum(I_zi);

disp('-------------------------')
if lr == 'l'
  disp('Left leg section')
elseif lr == 'r'
  disp('Right leg section')
end
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',mass)
fprintf('Volume:   %1.4f m^3\n',sum(v))
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [sum(m),sum(v),xc,yc,zc, Ip_x,Ip_y,Ip_z];

%% Plot

O17 = O16+[0;0;-l];
% axes:
plot_coord(O16,'index','13''');
plot_coord(O17,'index','14''');
 
opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};
 
for ii = indf
  ph = -ii*l/N; % plate height
  plot_elliptic_plate(O16+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
end

%% sideways paraboloids

h = 0.59*r;
a = r;
b = r;
n = 10;

nu = linspace(0,2*pi,n); % row
u = linspace(0,h,n)';    % column

x = a*sqrt(u/h)*cos(nu);
y = b*sqrt(u/h)*sin(nu);
z = u*ones(1,n);

surf(O17(1)+z-a-h,O17(2)+y,O17(3)+x,'facealpha',0.1,'edgealpha',0.2)
surf(O17(1)-z+a+h,O17(2)+y,O17(3)+x,'facealpha',0.1,'edgealpha',0.2)


end
