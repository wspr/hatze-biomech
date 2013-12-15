function person = segment_shoulder(person,S)

P = person.segment(2).origin+person.segment(S).offset;
R = person.segment(S).Rglobal;
i_m = person.sex;

gamma_1 = person.density.shoulder_lateral(i_m);
gamma_2 = person.density.shoulder_medial(i_m);
gamma_T = person.density.shoulder_cutout(i_m);
gamma_s = person.density.humerous(i_m);

l_t = person.meas{1}.length;
at1 = person.segment(1).a(1);
at4 = person.segment(1).a(4);
at5 = person.segment(1).a(5);

d   = person.meas{S}.all(1);
b   = person.meas{S}.all(2)/2;
b1  = person.meas{S}.all(3)/2;
z_h = person.meas{S}.all(4);

bt1 = person.meas{1}.depths(1)/2;
bt4 = person.meas{1}.depths(4)/2;


%% Calcs

j1 = 0.35*l_t-z_h;
j2 = 0.20*l_t-z_h;
h1 = 0.68*at5-at1;
d_z = j2-1.5*b1;

beta  = asin( d_z/d );
alpha = atan( h1/j1 );

d_x = d*cos(beta);
h_x = d_x - h1;
h_z = j1 - 2.5*b1 - d_z;
gamma = atan( h_z/h_x );

c2  = -(tan(beta)+tan(gamma));
c4  = -(b-1.42*b1)/h_x;
c1  = 2.5*b1  - d_x*c2;
c3  = 1.42*b1 - d_x*c4;
c5  = 0; % 0.35*l_t-z_h-h1*cot(alpha) = j1-h1/tan(alpha) = 0 by definition!
c6  = j1/h1 - d_z/d_x; % cot(alpha) - tan(beta)
c8  = bt4 - 0.15*l_t*(bt4-bt1)/j1;
c9  = (bt4-bt1)/j1;
c10 = (0.42*at4-at1)/j1;
c11 = (2*at1-1.42*at4)/(j1^2);

A1 = @(z) c1 + c2*z; % h1..d_x
B1 = @(z) c3 + c4*z;
% A2 = @(z) c5 + c6*z; % 0..h1
% B2 = @(z) b*sqrt(z/h1);

Oshoulder = P + person.segment(1).Rglobal*[0;0;-z_h-d_z-1.5*b1];
Oarm = Oshoulder + R*[ 0 ; 0; at1+d_x ]; % ??
person.segment(S).origin = Oshoulder;
person.segment(S+1).origin = Oarm;


% Mass 

v1 = 4/3*( ...
    c1*c3*h_x + (c2*c3 + c1*c4)*( (h1+h_x)^2 - h1^2 )/2 + ...
    c2*c4*( (h1+h_x)^3 - h1^3 )/3 ...
  );
v2 = 8/3*b*(1/3*c5*h1+1/5*c6*h1^2);
v_s = 2*pi*(b1/2)^3/3;

at = @(e) at4+c10*(e+0.15*l_t)+c11*(e+0.15*l_t).^2;
bt = @(e) c8-c9*e;
u  = @(e) (at1+(j2-e)*tan(alpha))/at(e);
fun = @(e) bt(e).*at(e).*(pi/2-u(e)*sqrt(1-u(e)^2)-asin(u(e)));
v_T = integral(fun,-0.15*l_t,j2);

m1  = gamma_1*v1;
m2  = gamma_2*v2;
m_s = gamma_s*v_s;
m_T = gamma_T*v_T;

v = v1+v2-v_s-v_T;
m = m1+m2-m_s-m_T;

% Mass centroid:
xc = 0;
yc = 0;

z2 = gamma_2*8/3*b*(1/5*c5*h1^2+1/7*c6*h1^3)/m2;
fun2 = @(e) bt(e).*at(e).*(pi/2-u(e).*sqrt(1-u(e)^2)-asin(u(e))).*e;
e_barm = ...
  m2*j2*(1-z2/h1) + ...
  -gamma_T*integral(fun2,-0.15*l_t,j2);
e_bar = e_barm/m;

fun3 = @(e) bt(e).*at(e).*(at1*(u(e).*(1-u(e)^2)^0.5-asin(u(e))-pi/2)+2*at(e).*(1-u(e)^2)^1.5/3);
JT = integral(fun3,j2,-0.15*l_t);
z_barm = 1.333*gamma_1*(c1*c3*((h1+h_x)^2-h1^2)/2+(c3*c2+c1*c4)*((h1+h_x)^3-h1^3)/3+c2*c4*((h1+h_x)^4-h1^4)/4)+8*gamma_2*b*(c5*h1^2/5+c6*h1^3/7)/3-m_s*(d_x-3*b1/16)-JT;
z_bar = z_barm/m;
theta7 = atan(e_bar/(d_x-z_bar));
zc = (at1+z_bar)/cos(theta7);

% principal moments of inertia; 
f1 = @(e)at(e)-at1-(j2-e)*tan(alpha);
f2 = @(e)((j2-e)*tan(alpha)+at1)/at(e);
f3 = c1*c3^3*h_x+(c2*c3^3+3*c1*c3^2*c4)*((h1+h_x)^2-h1^2)/2+3*(c1*c3*c4^2+c2*c3^2*c4)*((h1+h_x)^3-h1^3)/3+(c1*c4^3+3*c2*c3*c4^2)*((h1+h_x)^4-h1^4)/4+c2*c4^3*((h1+h_x)^5-h1^5)/5;
f4 = c3*c1^3*h_x+(c4*c1^3+3*c3*c1^2*c2)*((h1+h_x)^2-h1^2)/2+3*(c3*c1*c2^2+c4*c1^2*c2)*((h1+h_x)^3-h1^3)/3+(c3*c2^3+3*c4*c1*c2^2)*((h1+h_x)^4-h1^4)/4+c4*c2^3*((h1+h_x)^5-h1^5)/5;

% NB ambiguous: should c3 be swapped with c1, and c2 swapped with c4 (as
% above), or should c1 <-> c2 and c3 <-> c4 ? Hatze writes "f4 is the
% function f3 with c1 and c3, and c2 and c4 interchanged".

I_M = 4/3*gamma*(c1*c3*((h1+h_x)^3-h1^3)/3+(c2*c3+c1*c4)*((h1+h_x)^4-h1^4)/4+c2*c4*((h1-h_x)^5-h1^5)/5)-2*z_bar*1.333*gamma*(c1*c3*((h1+h_x)^2-h1^2)/2+(c2*c3+c1*c4)*((h1+h_x)^3-h1^3)/3+c2*c4*((h1+h_x)^4-h1^4)/4)+m1*z_bar^2;
I_s = m_s*(0.4*(b1/2)^2+(d_x-0.375*(b1/2)-z_bar)^2);
fun4 = @(e) 4*f1(e).*bt(e).*(1-f2(e).^2).^1.5/15+16*f1(e).^3.*bt(e).*(1-f2(e).^2)^0.5/175+4*f1(e).*bt(e).*(1-f2(e).^2).^0.5.*(z_bar-0.4*(at(e)-at1)-0.6*(j2-e)*tan(alpha)).^2/3;
fun5 = @(e) 16*f1(e).^3.*bt(e).*(1-f2(e).^2).^0.5/175+4*f1(e).*bt(e).*(1-f2(e).^2).^0.5.*(((z_h-0.4*(at(e)-at1)-0.6*(j2-e)*tan(alpha)).^2+(e-e_bar).^2)/3);
I_eT = integral (fun4,j2,-0.15*l_t);
I_nT = integral (fun5,j2,-0.15*l_t);

Ip_x = 0.2667*gamma*f3 + ...
  I_M - I_s + ...
  0.533*gamma*b^3*(c5*h1/5+c6*h1^2/7) + ...
  2.667*gamma*b*((c5*h1^3/7+c6*h1^4/9) + ...
  -2*z_bar*(c5*h1^2/5+c6*h1^3/7) + ...
  z_bar^2*(c5*h1/3+c6*h1^2/5)) + ...
  -I_eT;

Ip_y = 0.0914*gamma*f4+I_M+m1*e_bar^2-I_s-m_s*e_bar^2+0.1828*gamma*b*(c5^3*h1+3*c5^2*c6*h1^2/5+3*c5*c6^2*h1^3/7+c6^3*h1^4/9)+2.667*gamma*b*((1+(j2)^2/(h1^2))*(c5*h1^3/7+c6*h1^4/9)-2*(z_bar+((j2)^2-e_bar*(j2))/h1)*(c5*h1^2/5+c6*h1^3/7)+(z_bar^2+(e_bar-(j2))^2)*(c5*h1/3+c6*h1^2/5))-I_nT;
Ip_z = 0; % by definition

O1O7 = 0.8*l_t+e_bar*(d_x+at1)/(d_x-z_bar);
O7O8 = (d_x+at1)/cos(theta7);



%% finish

mass = m;
volume = v;

person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];

a10 = A1(d_x);
b10 = B1(d_x);
a1h = A1(h1);
b1h = B1(h1);

if person.plot || person.segment(S).plot

  if S == 7
    rcorr = [0 180 0];
    lr_sign = 1;
  else
    rcorr = [0 0 0];
    lr_sign = -1;
  end

  plot_parabolic_wedge(...
    Oarm,...
    [a10 b10],[a1h b1h],lr_sign*h_x,'skew',-h_z,'drop',-b1,...
    'rotate',R*rotation_matrix_zyx(rcorr),...
    'colour',person.segment(S).colour,...
    'face',[true false],...
    'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2))

  plot_parabolic_wedge(...
    Oarm+R*[0;0;-h_x],...
    [a1h b1h],0.001*[a1h b1h],lr_sign*h1,'skew',j1,'drop',-b1-h_z,...
    'face',[false true],...
    'rotate',R*rotation_matrix_zyx(rcorr),'colour',person.segment(S).colour,...
    'opacity',person.segment(S).opacity(1),'edgeopacity',person.segment(S).opacity(2))

end

end
