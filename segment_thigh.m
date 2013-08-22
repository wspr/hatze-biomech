function person = thigh(person,S,lr,h)

P = person.origin{S} + person.offset{S};
i_m = person.sex;

%% Thigh

N = 10; 
indf = 1:N;

%% Densities

nu = person.nu;
gamma = cellfun( @(x) x(i_m,nu), person.density.thigh );
gamma_0 = person.density.thigh_head(i_m,nu); 

%% Measurements

a = person.meas{S}.diam/2;
u = person.meas{S}.perim;
b = sqrt(((u/pi).^2)/2-a.^2);
l_1 = person.meas{S}.length;
h
h = person.meas{S}.length_long - l_1

person.meas{S}.a = a;
person.meas{S}.b = b;

%% Calculations

% Mass 
v = pi*a.*b*l_1/N;
m = gamma.*v;
v_0 = 2*pi*a(1)*b(1)*h/3;
m_0 = gamma_0*v_0;

volume = sum(v)+v_0;
mass   = sum(m)+m_0;

% Mass centroid:
xc = 0;
yc = 0;
zc = -(m_0*0.4*h+sum(h+l_1*(2*indf-1)/2/N))/mass;

% Moments of inertia:
I_x0 = m_0*(b(1)^2/4+0.0686*h^2);
I_y0 = m_0*(0.15*a(1)^2+0.0686*h^2);
I_z0 = m_0*(0.15*a(1)^2+b(1)^2/4);
I_xi = m.*(3*b.^2+(l_1/10)^2)/12;
I_yi = m.*(3*a.^2+(l_1/10)^2)/12;
I_zi = m.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = I_x0+m_0*(-0.4*h-zc).^2+sum(I_xi+m.*(h+l_1*(2*indf-1)/20+zc).^2);
Ip_y = I_y0+m_0*(-0.4*h-zc).^2+sum(I_yi+m.*(h+l_1*(2*indf-1)/20+zc).^2);
Ip_z = I_z0+sum(I_zi);

Q = P+[0;0;-l_1-h];
person.origin{S+1} = Q;

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot

opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
 
for ii = indf
  ph = l_1-ii*l_1/N; % plate height
  plot_elliptic_plate(Q+[0;0;ph],[a(ii) b(ii)],l_1/N,opt{:})
end

plot_hoof(P-[0;0;h],a(1),b(1),h,opt{:})

end

end
