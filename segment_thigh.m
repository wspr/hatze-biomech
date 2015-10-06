function person = segment_thigh(person,S)

P = person.segment(S).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;
N = person.segment(S).Ncalc;
PI = person.const.pi;

ind = 1:N;

L = person.meas{S}.length;
a = person.resample(person,S,person.meas{S}.diam)/2;
u = person.resample(person,S,person.meas{S}.perim);
b = person.solve_ellipse(a,u);

gamma = person.resample(person,S,cellfun( @(x) x(person.sex,person.nu), person.density.thigh ));
gamma_0 = person.density.thigh_head(person.sex,person.nu); 

% the height of the hoof is calculated per-thigh; for the pelvic region,
% however, it is set to a fixed height! Therefore, the hoof does not
% fit perfectly inside of its cutout in the pelvis; it is unclear whether
% this is intentional or not (in Hatze's eyes).
%h = person.meas{S}.length_long - L;

person.meas{S}.a = a;
person.meas{S}.b = b;

%% Hoof stuff from abdomino-pelvic section
h_hoof_l = person.meas{12}.length_long - person.meas{12}.length;
h_hoof_r = person.meas{15}.length_long - person.meas{15}.length;

dmean = (h_hoof_l + h_hoof_r)/2;

if abs((h_hoof_l-h_hoof_r)/dmean) > 0.07
   if h_hoof_r > h_hoof_l
       h_hoof_r = h_hoof_l + 0.07*dmean;
   elseif h_hoof_r < h_hoof_l
       h_hoof_l = h_hoof_r + 0.07*dmean;
   else
       h_hoof_l=h_hoof_r;    
   end
end

if S==12
   h=h_hoof_l;
else 
   h=h_hoof_r;
end

%% Calculations

% Mass
v = PI*a.*b*L/N;
m = gamma.*v;
v_0 = 2*PI*a(1)*b(1)*h/3; 
m_0 = gamma_0*v_0;

volume = sum(v)+v_0;
mass   = sum(m)+m_0;

% Mass centroid:
xc = 0;
yc = 0;
zc = -(m_0*0.4*h+sum(m(ind).*(h+L*(2*ind-1)/20)))/mass;

% Moments of inertia:
I_x0 = m_0*(1/4*b(1)^2+0.0686*h^2);
I_y0 = m_0*(0.15*a(1)^2+0.0686*h^2);
I_z0 = m_0*(0.15*a(1)^2+1/4*b(1)^2);

I_xi = m.*(3*b.^2+(L/N)^2)/12;
I_yi = m.*(3*a.^2+(L/N)^2)/12;
I_zi = m.*(a.^2+b.^2)/4;

% Principal moments of inertia (w.r.t centroid);
Ip_x = I_x0 + m_0*(-0.4*h-zc).^2 + sum(I_xi+m.*(h+L*(2*ind-1)/20+zc).^2);
Ip_y = I_y0 + m_0*(-0.4*h-zc).^2 + sum(I_yi+m.*(h+L*(2*ind-1)/20+zc).^2);
Ip_z = I_z0 + sum(I_zi);

%Principal moments of inertia w.r.t local systems origin
PIOX = Ip_x+mass*zc^2;
PIOY = Ip_y+mass*zc^2;
PIOZ = Ip_z;

%Coordinates of origin of axes - needs to be updated with new
%abdomino-pelvic section
%OX12 = -(person.meas{11}.all(8)-person.meas{S}.all(1))
% DT=person.meas{11}.r-person.meas{11}.all(18);
%OY12 = DT*cos(person.segment(11).theta)+(h-person.meas{11}.length)*sin(person.segment(11).theta)
%OZ12 = (h-person.meas{11}.length)*cos(person.segment(11).theta)-DT*sin(person.segment(11).theta)

Q = P+person.segment(S).Rglobal*[0;0;-L-h];
person.segment(S+1).origin = Q;

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot

  opt  = {'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2),'colour',person.segment(S).colour};

  for ii = ind
    ph = L-ii*L/N; % plate height
    plot_elliptic_plate(Q+person.segment(S).Rglobal*[0;0;ph],[a(ii) b(ii)],L/N,opt{:},'rotate',R)
  end

  plot_hoof(P-R*[0;0;h],a(1),b(1),h,opt{:},'rotate',R)

end

end
