function person = segment_thigh(person,S)

P = person.origin{S} + person.offset{S};
person.segment(S).Rglobal = person.segment(11).Rglobal*person.segment(S).Rlocal;
N = person.segment(S).Ncalc;

i_m = person.sex;
ind = 1:N;

l_1 = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2; 
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

gamma = person.resample(person,S,cellfun( @(x) x(person.sex,person.nu), person.density.thigh ));
gamma_0 = person.density.thigh_head(i_m,person.nu); 

% the height of the hoof is calculated per-thigh; for the pelvic region,
% however, it is set to a fixed height! Therefore, the hoof does not
% fit perfectly inside of its cutout in the pelvis; it is unclear whether
% this is intentional or not (in Hatze's eyes).
h = person.meas{S}.length_long - l_1;

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
zc = -(m_0*0.4*h+sum(h+l_1*(ind-1/2)/N))/mass;

% Moments of inertia:
I_x0 = m_0*(b(1)^2/4+0.0686*h^2);
I_y0 = m_0*(0.15*a(1)^2+0.0686*h^2);
I_z0 = m_0*(0.15*a(1)^2+b(1)^2/4);
I_xi = m.*(3*b.^2+(l_1/N)^2)/12;
I_yi = m.*(3*a.^2+(l_1/N)^2)/12;
I_zi = m.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = I_x0+m_0*(-0.4*h-zc).^2+sum(I_xi+m.*(h+l_1*(ind-1/2)/N+zc).^2);
Ip_y = I_y0+m_0*(-0.4*h-zc).^2+sum(I_yi+m.*(h+l_1*(ind-1/2)/N+zc).^2);
Ip_z = I_z0+sum(I_zi);

Q = P+person.segment(S).Rglobal*[0;0;-l_1-h];
person.origin{S+1} = Q;

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot
  
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = ind
    ph = l_1-ii*l_1/N; % plate height
    plot_elliptic_plate(Q+[0;0;ph],[a(ii) b(ii)],l_1/N,opt{:})
  end
  
  plot_hoof(P-[0;0;h],a(1),b(1),h,opt{:})
  
end

end
