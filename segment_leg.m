function person = segment_leg(person,S)

P = person.origin{S} + person.offset{S};
i_m = person.sex;

%% Leg

N = 10; 
ind = 1:N;

%% Densities

nu = person.nu; % subcutaneous fat indicator;
gamma = cellfun( @(x) x(i_m), person.density.leg );
gamma_b = person.density.ankle;

%% Measurements

a = person.meas{16}.diam/2; 
u = person.meas{16}.perim;
b = sqrt (((u(ind)/pi).^2)/2-a(ind).^2); 
l = person.meas{16}.length;
r = person.meas{16}.ankle;
r_a = 0.59*r;

%% Calculations

% Volume
v      = pi*a.*b*l/N;
v_b    = pi*r_a/2*r^2;
volume = sum(v) + 2*v_b;

% Mass
m    = gamma.*v;
m_b  = gamma_b*v_b; % ankle
mass = sum(m) + 2*m_b;

% Mass centroid:
xc = 0;
yc = 0;
zc = -2*m_b*l - sum(m.*l.*(2*ind-1)/20)./mass;

% Moments of inertia:
I_xi = m.*(3*b.^2+(l/10)^2)/12;
I_yi = m.*(3*a.^2+(l/10)^2)/12;
I_zi = m.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = 2*m_b*(0.33*r^2+(l+zc).^2)+sum(I_xi+m.*(l*(2*ind-1)/20+zc).^2);
Ip_y = 2*m_b*(0.1859*r^2+(l+zc).^2+(a(10)+0.196*r^2)) + sum(I_xi+m.*(l*(2*ind-1)/20+zc).^2);
Ip_z = 2*m_b*(0.1859*r^2+(a(10)+0.196*r)^2)+ sum(I_zi);

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

Q = P+[0;0;-l];
person.origin{S+1} = Q;

%% Plot

if person.plot
  
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = ind
    ph = -ii*l/N; % plate height
    plot_elliptic_plate(P+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
  end
  
  %% sideways paraboloids
  
  a = r;
  b = r;
  n = 10;
  
  nu = linspace(0,2*pi,n); % row
  u = linspace(0,r_a,n)';    % column
  
  x = a*sqrt(u/r_a)*cos(nu);
  y = b*sqrt(u/r_a)*sin(nu);
  z = u*ones(1,n);
  
  surf(Q(1)+z-a-r_a,Q(2)+y,Q(3)+x,'facealpha',person.opacity{S}(1),'edgealpha',person.opacity{S}(2),'facecolor',person.color{S})
  surf(Q(1)-z+a+r_a,Q(2)+y,Q(3)+x,'facealpha',person.opacity{S}(1),'edgealpha',person.opacity{S}(2),'facecolor',person.color{S})
  
end

end
