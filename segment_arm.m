function person = segment_arm(person,S)

P = person.origin{S} + person.offset{S};
i_m = person.sex;

%%  Arm

N = person.segment(S).N; 
indu = 1:N;
 
%% Densities

gamma_0 = person.density.humerous(i_m); 
gamma = cellfun( @(x) x(i_m), person.density.arm );

%% Measurements

a = person.meas{S}.diam/2;
u = person.meas{S}.perim;
b = person.solve_ellipse(a,u);
l = person.meas{S}.length;

Q = P+[0;0;-l];
person.origin{S+1} = Q;

%% Calculations

% volume
v      = pi*a.*b*l/N; 
v_0    = 2*pi*((b(1)/2)^3)/3;
volume = sum(v) + v_0;

% mass
m    = gamma.*v;
m_0  = gamma_0*v_0;
mass = m_0 + sum(m); 
 
% Mass centroid:
xc = 0;
yc = 0;
zc = (m_0*0.375*(b(1)/2) - sum(m.*(2*indu-1).*l)/2/N)/mass;
 
% Moments of inertia:
I_x = mass.*(3*(b.^2)+(l/10).^2)/12; 
I_y = mass.*(3*(a.^2)+(l/10).^2)/12;
I_z = mass.*(3*(b.^2)+(b/10).^2)/12;
 
% principal moments of inertia; 
Ip_x = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_x) +sum(m.*(l*(2*indu-1)/20+zc).^2);
Ip_y = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_y) +sum(m.*(l*(2*indu-1)/20+zc).^2);
Ip_z= m_0*((b(1)./2)^2)/5 + sum(I_z);
 
person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot,...

if person.plot
  
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = indu
    ph = l-ii*l/N; % plate height
    plot_elliptic_plate(Q+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
  end
  
  % the hemishpere
  plot_sphere(P, b(1)/2, 'longrange',[0 1],opt{:})
  
end

end