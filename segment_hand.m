 function person = segment_hand(person,S)

if S == 6
  lr = -1;
else
  lr = 1;
end

P = person.segment(S).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;

PI = person.const.pi;

r = person.meas{S}.all(1);
h = person.meas{S}.all(2);

a10 = person.meas{S-1}.a(end);
b10 = person.meas{S-1}.b(end);

hh = 3.4/9*h;

density = person.density.hand;

%% Calculations
% Volume and mass

v_p = b10*hh*(2.172*a10+1.172*h);
%v_c = PI*h^2*(2*r-h/4)/8;  % Hatze 79
v_c = h^2*(0.7854*r-0.098*h);   %fortran code (decimal values of the above equation)
v_T = (h^3)/16;

m_p = density*v_p;
m_c = density*v_c;
m_T = density*v_T;

volume = v_p + v_c + v_T;
mass   = m_p + m_c + m_T;

% Mass centroid  (w.r.t original segment axes)

xc1= 0.424*(r^3-(r-h/4)^3)/(r^2-(r-h/4)^2);  %fortran code has 0.424, Hatze 79 has 0.4244
zp1=(density*b10*(hh^2)*(0.753*a10+0.797*h))/m_p;  %fortran code has 0.753 instead of 0.7533 and 0.797 instead of 0.7967

xc = lr*(m_c*xc1-m_T*(r-h/4))/mass;
yc = (m_T*3*h/8)/mass;
zc = -(m_p*zp1+m_c*(r+hh)+m_T*(r/3+hh))/mass;

% Moments of inertia

Ig_xc=m_c*(r^2+(r-h/4)^2+(h^2)/3)/4;
Ig_yc=m_c*((r^2+(r-h/4)^2)/2-xc1^2);
Ig_zc=Ig_xc-m_c*xc1^2;

c1 = 0.26/hh;
c2 = (h/2-a10)/hh;

I_x1=1.333*density*b10*((a10^3*hh)+(3*a10^2*c2+c1*a10^3)*(hh^2)/2+...
    3*(a10*c2^2+c1*a10^2*c2)*(hh^3)/3+...
    (c2^3+3*c1*a10*c2^2)*(hh^4)/4+...
    c1*c2^2*((hh^5)/5));
I_y2=1.333*density*b10^3*(a10*hh+(3*a10*c1+c2)*(hh)^2/2+...
    3*(a10*c1^2*+c1*c2)*(hh)^3/3+...
    (a10*c1^3+3*c2*c1^2)*(hh)^4/4+...
    (c1^3*c2*(0.378*h)^5)/5);
Ip=density*b10*(hh^3)*(0.385*a10+0.6040*h)-m_p*zp1^2;

I_x = I_x1+Ip+m_p*(yc^2+(zp1+zc)^2)+Ig_xc+...
    m_c*(yc^2+(r+hh+zc)^2)+...
    m_T*(((3*h)/8-yc)^2+(r/3+hh+zc)^2);
I_z = I_x1+I_y2+m_p*(xc^2+yc^2)+Ig_xc+...
    m_c*(yc^2+(xc1-abs(xc))^2-xc1^2)+...
    m_T*((abs(xc)+r-h/4)^2+(3*h/8-yc)^2);
I_zx = lr*density*h*((r^4-(r-h/4)^4)*(sin(h/(r-h/8)))^2/32-((r^3-(r-h/4)^3)/3)*...
    ((xc/4)*sin(h/(r-h/8))+(zc+r+hh)*(2.25-0.25*cos(h/(r-h/8))))+...
    0.5*(r^2-(r-h/4)^2)*xc*(zc+r+hh)*(PI+0.25*h/(r-h/8)))+lr*m_p*(-xc)*(-zc-zp1);

% principal moments of inertia;
Ip_y = I_y2+Ip+m_p*(xc^2+(zp1+zc)^2)+Ig_yc+m_c*((xc1-abs(xc))^2+(r+hh+zc)^2)+m_T*((abs(xc)+r-h/4)^2+(r/3+hh+zc)^2); %Note xc^2
Ip_x =(I_z+I_x)/2-sqrt(1/4*(I_z-I_x)^2+I_zx^2);
Ip_z =(I_z+I_x)/2+sqrt(1/4*(I_z-I_x)^2+I_zx^2);

theta = atan(I_zx/(I_x-Ip_z));

%centroid w.r.t local coordinate systems (since principal axes differ from
%original segment axes)
xbc=xc*cos(theta)-zc*sin(theta);
ybc=yc;
zbc=zc*cos(theta)+xc*sin(theta);

%prinicpal moments of inertia w.r.t local systems origin
PIOX=Ip_x+mass*(ybc^2+zbc^2);
PIOY=Ip_y+mass*(xbc^2+zbc^2);
PIOZ=Ip_z+mass*(xbc^2+ybc^2);

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc]; 
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z]; 
person.segment(S).theta = theta;  

%%
if person.plot || person.segment(S).plot

  opt = {'colour',person.segment(S).colour,'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2)};

  plot_rect_prism(P,[2*a10 2*b10],[h 2.52*b10],-hh,'rotate',R,opt{:});

  if S == 6 % left
    s = h;
    t = 0;
    RR = R*person.cardan_rotation([0 0 180]);
  else
    s = 0;
    t = 3*h/4;
    RR = R;
  end

  plot_cylinder_hollow(P+R*[0;s-h/2;-hh-r],r-h/4,r,h,[0 180],'N',10,'rotate',RR,opt{:})
  plot_cylinder_hollow(P+R*[0;s+t-h/2;-hh-r],r-h/4,r,h/4,[0 -180/PI*h/r],'N',10,'rotate',RR,opt{:})

end

end
