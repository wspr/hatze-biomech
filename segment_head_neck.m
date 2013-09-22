function person = segment_head_neck(person,S)

%% Head_neck

P = person.origin{S}+person.offset{S};

i_m = person.sex;

head_width  = person.meas{2}.all(1);
head_depth  = person.meas{2}.all(2);
head_height = person.meas{2}.all(3);
neck_height = person.meas{2}.all(4);

%% Densities:

gamma_e = person.density.head(i_m);
gamma_c = person.density.neck(i_m);

%% Measurements

a = head_width/2;
b = 1.08*head_depth/2;
c = 1.04*head_height/2;
a_1 = person.meas{1}.all(1)/2;
b_1 = person.meas{1}.all(9)/2;
h = neck_height;
k = 0.30*c;

%% Calculations

% Volume 
v_e = 4.66493*a*b*c;
v_c = pi*a_1*b_1*h;
v_p = 4.5708*b*c*a_1*(1.5708-b_1*(c-k)/(c*b)-asin((c-k)/c));
volume = v_e + v_c + v_p;

% Mass 
m_e = gamma_e*v_e;
m_c = gamma_c*v_c;
m_p = gamma_e*v_p;
mass = m_e + m_c + m_p;

% Mass centroid:
xc = 0;
yc = 0;
zc = (m_e*(c+h-k)+m_c*h/2-m_p*(h-k/3))/mass;

% Moments of inertia:
I_xe = mass*(0.19473*b^2+0.23511*c^2);
I_ye = mass*(0.211*a^2+0.23511*c^2);
I_ze = mass*(0.211*b^2+0.19473*b^2);
I_xc = mass*(3*b^2+h^2)/12;
I_yc = mass*(3*a^2+h^2)/12;
I_zc = mass*(a^2+b^2)/4;

% principal moments of inertia; 
Ip_x = I_xe+m_e*(c+h-k-zc)^2+I_xc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_y = I_ye+m_e*(c+h-k-zc)^2+I_yc+m_c*(zc-h/2)^2-m_p*(zc-h+k/3)^2;
Ip_z = I_ze+I_zc-m_p*(a_1^2+b_1^2)/6;

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

%% Plot

if person.plot || person.segment(S).plot
  
  plot_elliptic_plate(P,[a_1 b_1],h,'colour',person.color{S},...
    'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2));
  
  zf = @(x,y) c*sqrt(1-(y/b).^2).*(1-(x./a).^8);
  
  N = 21;
  t = linspace(-pi,pi,N);
  d = linspace(0,1,N);
  
  [tt,dd] = meshgrid(t,d);
  
  ho = h + c/2 + k;
  
  x = dd*a.*cos(tt);
  y = dd*b.*sin(tt);
  z = zf(x,y);
  
  opt = {'facecolor',person.color{S},...
    'facealpha',person.opacity{S}(1),...
    'edgealpha',person.opacity{S}(2)};
  
  surf(P(1)+x,P(2)+y,P(3)+ho+z,opt{:})
  surf(P(1)+x,P(2)+y,P(3) + ho-z,opt{:})
  
  x = a*sin(t);
  y = b*cos(t);
  z = zf(x,y);
  surf(P(1)+[x;x],P(2)+[y;y],P(3)+ho+[z;-z],opt{:})
  
end

end
