function person = segment_arm(person,S)

P = person.origin{S} + person.offset{S};
R = person.segment(S).Rglobal;

N = person.segment(S).Ncalc;
ind = 1:N;

L = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2; 
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

gamma = person.resample(person,S,cellfun( @(x) x(person.sex), person.density.arm ));
gamma_0 = person.density.humerous(person.sex); 

Q = P+person.segment(S).Rglobal*[0;0;-L];
person.origin{S+1} = Q;

%% Calculations

% volume
v      = pi*a.*b*L/N; 
v_0    = 2*pi*((b(1)/2)^3)/3;
volume = sum(v) + v_0;

% mass
m    = gamma.*v;
m_0  = gamma_0*v_0;
mass = m_0 + sum(m); 
 
% Mass centroid:
xc = 0;
yc = 0;
zc = (m_0*0.375*(b(1)/2) - sum(m.*(2*ind-1).*L)/2/N)/mass;
 
% Moments of inertia:
I_x = m.*(3*(b.^2)+(L/N).^2)/12; 
I_y = m.*(3*(a.^2)+(L/N).^2)/12;
I_z = m.*(a.^2+b.^2)/4;
 
% principal moments of inertia; 
Ip_x = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_x) +sum(m.*(L*(ind-1/2)/N+zc).^2);
Ip_y = m_0*((0.259*((b(1)/2)^2))+((0.375*b(1))/(2-zc))^2)+sum(I_y) +sum(m.*(L*(ind-1/2)/N+zc).^2);
Ip_z= m_0*((b(1)./2)^2)/5 + sum(I_z);
 
person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot...

if person.plot || person.segment(S).plot
  
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = ind
    ph = L-ii*L/N; % plate height
    plot_elliptic_plate(Q+R*[0;0;ph],[a(ii) b(ii)],L/N,opt{:},'rotate',R)
  end
  
  % the hemishpere
  plot_sphere(P, b(1)/2, 'longrange',[0 1],'rotate',R,opt{:})
  
end

end