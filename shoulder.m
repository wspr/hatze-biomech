function [O8,O7] = shoulder(O2,i_m,lr,a_thorax,l_t,d_arm)
%%

meas = [0.140 0.176 0.090 0.025];

b = meas(1)/2;
d = meas(2);
% meas(3) ?
z_h = meas(4);

b1 = d_arm(1)/4;

j1 = 0.35*l_t-z_h;
h1 = 0.68*a_thorax(5)-a_thorax(1);
d_z = 0.2*l_t-z_h-1.5*b1;

beta = asin( d_z/d );
alpha = atan( h1/j1 );

d_x = d*cos(beta);
hh = d_x - h1;
gamma = atan( (j1 - 2.5*b1 - d_z)/hh );

c2 = (tan(beta)+tan(gamma));
c4 = (b-1.42*b1)/hh;
c1 = 2.5*b1  + d_x*c2;
c3 = 1.42*b1 + d_x*c4;

A1 = @(z) c1 - c2*z; % h1..d_x
B1 = @(z) c3 - c4*z;

% c5 = j1-h1/tan(alpha); % these two terms are equal!! ( i think )
c5 = 0;
c6 = 1/tan(alpha)-tan(beta);

A2 = @(z) c5 + c6*z; % 0..h1
B2 = @(z) b*sqrt(z/h1);

O7 = O2 + [0;0;-z_h-d_z-1.25*b1];

if lr == 'r'
  lr_sign = 1;
elseif lr == 'l'
  lr_sign = -1;
end

O8 = O7 + [ lr_sign*(a_thorax(1)+d_x) ; 0; 0]

start = O2+[lr_sign*a_thorax(1); 0; -z_h];
shtop = start+[lr_sign*d_x;0;-d_z];

a10 = A1(d_x);
b10 = B1(d_x);
a1h = A1(h1);
b1h = B1(h1);

a20 = A2(0);
b20 = B2(0);
a2h = A2(h1);
b2h = B2(h1);

shift = -hh*sin(gamma);
plot_parabolic_wedge([O8(1);0;O2(3)-z_h-d_z-a10],[a10 b10],[a1h b1h],lr_sign*hh,'t',shift,'rotate',[0 -90 0])
plot_parabolic_wedge([O8(1)-lr_sign*hh;0;O2(3)-z_h-d_z-a10+shift],[a2h b2h],0.001*[a2h b2h],lr_sign*h1,'t',j1,'rotate',[0 -90 0])

X = h1/tan(alpha);
Y = hh*tan(gamma);

end
