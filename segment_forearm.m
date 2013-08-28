function person = segment_forearm(person,S)

P = person.origin{S} + person.offset{S}(:);
person.segment(S).Rglobal = person.segment(S-1).Rglobal*person.segment(S).Rlocal;
N = person.segment(S).Ncalc;
ind = 1:N;

l = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2; 
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

gamma = person.resample(person,S,cellfun( @(x) x(person.sex), person.density.forearm ));

person.meas{S}.a = a;
person.meas{S}.b = b;

Q = P+person.segment(S).Rglobal*[0;0;-l];
person.origin{S+1} = Q; 

% Mass 
v = pi*a.*b*l/N; % volume of each forearm disk
m = gamma.*v; % mass of each forearm disk

mass = sum(m);
volume = sum(v);

% Mass centroid:
xc = 0;
yc = 0;
zc = -sum(m.*(ind-1/2).*l)/(N*mass);

% Moments of inertia:
I_x = m.*(3*(b.^2)+(l/N).^2)/12; 
I_y = m.*(3*(a.^2)+(l/N).^2)/12;
I_z = m.*(a.^2+b.^2)/4;

% principal moments of inertia; 
Ip_x = sum(I_x) + sum(m.*(l*(ind-1/2)/N+zc).^2);
Ip_y = sum(I_y) + sum(m.*(l*(ind-1/2)/N+zc).^2);
Ip_z=sum(I_z);

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot
  
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = ind
    ph = l-ii*l/N; % plate height
    plot_elliptic_plate(Q+[0;0;ph],[a(ii) b(ii)],l/N,opt{:})
  end
  
end

end

