function person = segment_shoulder(person,S)

P = person.origin{2}+person.offset{S};
person.segment(S).Rglobal = person.segment(2).Rglobal*person.segment(S).Rlocal;
i_m = person.sex;

gamma_1 = person.density.shoulder_lateral(i_m);
gamma_2 = person.density.shoulder_medial(i_m);
gamma_T = person.density.shoulder_cutout(i_m);
gamma_s = person.density.humerous(i_m);

l_t = person.meas{1}.length;
at1 = person.meas{1}.widths(1)/2;
at5 = person.meas{1}.widths(5)/2;

b = person.meas{S}.all(2)/2;
d = person.meas{S}.all(1);
% meas(3) ?
z_h = person.meas{S}.all(4);

b1 = person.solve_ellipse(person.meas{S+1}.diam(1)/2,person.meas{S+1}.perim(1))/2;

j1 = 0.35*l_t-z_h;
h1 = 0.68*at5-at1;
d_z = 0.2*l_t-z_h-1.5*b1;

beta = asin( d_z/d );
alpha = atan( h1/j1 );

d_x = d*cos(beta);
hh = d_x - h1;
gamma = atan( (j1 - 2.5*b1 - d_z)/hh );

X = h1/tan(alpha);
Y = hh*tan(gamma);

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

O7 = P + person.segment(1).Rglobal*[0;0;-z_h-d_z-1.25*b1];

if S == 3 % left
  lr_sign = 1;
elseif S == 7 % right
  lr_sign = -1;
end

O8 = P + person.segment(1).Rglobal*[ lr_sign*(at1+d_x) ; 0; -z_h-d_z-1.25*b1];
person.origin{S} = O7;
person.origin{S+1} = O8;

start = P+[lr_sign*at1; 0; -z_h];
shtop = start+[lr_sign*d_x;0;-d_z];

a10 = A1(d_x);
b10 = B1(d_x);
a1h = A1(h1);
b1h = B1(h1);

a20 = A2(0);
b20 = B2(0);
a2h = A2(h1);
b2h = B2(h1);


if person.plot
  
  s = -hh*sin(gamma);
  plot_parabolic_wedge(...
    [O8(1);O8(2);P(3)-z_h-d_z-a10],...
    [a10 b10],[a1h b1h],lr_sign*hh,'t',s,...
    'rotate',person.segment(S).Rglobal,...
    'colour',person.color{S},...
    'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2))
  plot_parabolic_wedge(...
    [O8(1)-lr_sign*hh;O8(2);P(3)-z_h-d_z-a10+s],...
    [a2h b2h],0.001*[a2h b2h],lr_sign*h1,'t',j1,...
    'rotate',person.segment(S).Rglobal,'colour',person.color{S},...
    'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2))
  
end

end
