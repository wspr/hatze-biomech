function person = segment_hand(person,S)

if S == 6
  lr = -1;
else
  lr = 1;
end

P = person.segment(S).origin + person.segment(S).offset;
R = person.segment(S).Rglobal;

r = person.meas{S}.all(1);
h = person.meas{S}.all(2);

a10 = person.meas{S-1}.a(end);
b10 = person.meas{S-1}.b(end);

hh = 0.378*h;

density = person.density.hand;

%% Volume and mass

v_p = b10*(0.378*h)*(2.172*a10+1.172*h);
v_c = pi*h^2*(2*r-h/4)/8;
v_T = h^3/16;

m_p = density*v_p;
m_c = density*v_c;
m_T = density*v_T;

volume = v_p + v_c + v_T;
mass   = m_p + m_c + m_T;

%% Mass centroid

xc1=0.4244*(r^3-(r-h/4)^3)/(r^2-(r-h/4)^2);
zp1=r*b10*(0.378*h)^2*(0.7533*a10+0.7967*h)/m_p;

xc = lr*(m_c*xc1-m_T*(r-h/4))/mass;
yc = m_T*3*h/8/mass;
zc = -(m_p*zp1+m_c*(r+0.378*h)+m_T*(r/3+0.378*h))/mass;

%% Moments of inertia

Ig_xc=m_c*(r^2+(r-h/4)^2+h^2/3)/4;
Ig_yc=m_c*(r^2+(r-h/4)^2)/2-xc1^2;
Ig_zc=Ig_xc-m_c*xc1^2;

c1 = 0.26/0.378*h;
c2 = (h/2-a10)/(0.378*h);

I_x1=1.333*r*b10*(a10^3*0.378*h)+(3*a10^2*c2+...
    c1*a10^3)*(0.378*h)^2/2+3*(a10*c2^2+...
    c1*a10^2*c2)*(0.378*h)^3/3+(c2^3+...
    3*c1*a10*c2^2)*(0.378*h)^4/4+...
    c1*c2^2*(0.378*h)^5/5;
I_y2=1.333*r*b10^3*(a10*(0.378*h)+...
    (3*a10*c1+c2*(0.378*h)^2/2+...
    3*(a10*c1^2*+c1*c2)*(0.378*h)^3/3+...
    (a10*c1^3+3*c1*c2^2)*(0.378*h)^4/4+...
    c1^3*c2*(0.378*h)^5/5));
Ip=r*b10*(0.378*h)^3*(0.385*a10+0.604*h)-m_p*zp1^2;

I_x = I_x1+Ip+m_p*(yc^2+(zp1+zc^2))+Ig_xc+m_c*(yc^2+(r+0.378*h+zc)^2)+m_T*((3*h/8-yc)^2+(r/3+0.378*h+zc)^2);
I_z=I_x1+I_y2+m_p*(xc^2+yc^2)+Ig_xc+m_c*(yc^2+(xc1-abs(xc))^2-xc1^2)+m_T*((abs(xc)+r-h/4)^4+(3*h/8-yc)^2);
I_xz= lr*r*h*((r^4-(r-h/4)^4))*(sin(h/(h/8-r)))^2/32-((r^3-(r-h/4)^3)*3)*((xc/4)*sin(h/(h/8-r)))+...
    (zc+r+0.378*h)*(2.25-0.25*cos(h/(h/8-r)))+...
    0.5*(r^2-(r-h/4)^2*xc*(zc+r+0.378*h)*(pi+0.25*h/(h/8-r)))+lr*m_p*(-xc)*(-zc-zp1);

% principal moments of inertia;
Ip_y =I_y2+m_p*(xc^2+(zp1+zc))^2+Ig_yc+m_c*(xc1-abs(xc))^2+(r+0.378*h+zc)^2+m_T*((abs(xc)+r-h/4)^2+(r/3+0.378*h+zc)^2);
Ip_x =(I_x+I_z)/2+sqrt(1/4*(I_x-I_z)^2+I_xz^2);
Ip_z =(I_x+I_z)/2-sqrt(1/4*(I_x-I_z)^2+I_xz^2);

theta = atan(I_xz/(I_z-Ip_x)); 

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
  plot_cylinder_hollow(P+R*[0;s+t-h/2;-hh-r],r-h/4,r,h/4,[0 -180/pi*h/r],'N',10,'rotate',RR,opt{:})

end

end
