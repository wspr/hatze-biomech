function [calcs, O6] = forearm(O5,i_m,lr,forearm_diameters,forearm_perimeters,forearm_length)

%% Forearm

N = 10; 
indf = 1:N;

%% Densities

gamma_i1 = @(ii,i_m) (1160-60*ii)*(1+0.0213*i_m); % for i = 1,2;
gamma_i2 = @(ii,i_m) (1034.2-2.86*ii)*(1+0.0213*i_m); % for i = 3,4,5,6,7,8;
gamma_i3 = @(i_m) 1204.29*(1+0.0213*i_m); % for i =9,10;
 
gamma_i = @(i_m) [gamma_i1([1 2],i_m) gamma_i2(3:8,i_m) gamma_i3(i_m) gamma_i3(i_m)];

%% Calculations

% 1:10 diameters 2ai 
a = forearm_diameters/2; % ai 
u = forearm_perimeters;
l = forearm_length;

% 11:20 perimieter ui
b =sqrt (((u(indf)/pi).^2)/2-a(indf).^2); % b

% Mass 
v = pi*a(indf).*b(indf)*l/N; % volume of each forearm disk
m = gamma_i(i_m).*v; % mass of each forearm disk

% Mass centroid:
xc = O5(1);
yc = O5(1);
zc = O5(3) - sum(sum(m(indf).*(2*indf-1).*l)/20*m);

% Moments of inertia:
I_x = m.*(3*(b.^2)+(l/10).^2)/12; 
I_y = m.*(3*(a.^2)+(l/10).^2)/12;
I_z = m.*(3*(b.^2)+(b/10).^2)/12;

% principal moments of inertia; 
Ip_x = sum(I_x) +sum(m.*(l*(2*indf-1)/20+zc).^2);
Ip_y = sum(I_y) +sum(m.*(l*(2*indf-1)/20+zc).^2);
Ip_z=sum(I_z);

mass = sum(m);
volume = sum(v);

disp('-------------------------')
if lr == 'l'
disp('Left forearm section')
elseif lr == 'r'
disp('Right forearm section')
end
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',mass)
fprintf('Volume:   %1.4f m^3\n',volume)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [mass,volume,xc,yc,zc,Ip_x,Ip_y,Ip_z];

%% Plot

O6 = O5+[0;0;-l];
 
opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};
 
for ii = indf
  ph = l-ii*l/N; % plate height
  plot_elliptic_plate(O6+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
end

end

