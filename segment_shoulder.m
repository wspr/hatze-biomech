function person = segment_shoulder(person,S)

P = person.segment(2).origin+person.segment(S).offset;
R = person.segment(S).Rglobal;
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
% alpha = atan( h1/j1 );

d_x = d*cos(beta);
h_x = d_x - h1;
h_z = j1 - 2.5*b1 - d_z;
gamma = atan( h_z/h_x );

c2 = (tan(beta)+tan(gamma));
c4 = (b-1.42*b1)/h_x;
c1 = 2.5*b1  + d_x*c2;
c3 = 1.42*b1 + d_x*c4;

A1 = @(z) c1 - c2*z; % h1..d_x
B1 = @(z) c3 - c4*z;

% Unnecessary for plotting:
%
% c5 = j1-h1/tan(alpha); % these two terms are equal!! hatze is silly
% c5 = 0;
% c6 = 1/tan(alpha)-tan(beta);
% A2 = @(z) c5 + c6*z; % 0..h1
% B2 = @(z) b*sqrt(z/h1);

Oshoulder = P + person.segment(1).Rglobal*[0;0;-z_h-d_z-1.5*b1];
Oarm = Oshoulder + R*[ 0 ; 0; at1+d_x ];
person.segment(S).origin = Oshoulder;
person.segment(S+1).origin = Oarm;

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
