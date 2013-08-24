
% not sure if I want to use these or not:
abthor         = 1;
headneck       = 2;
shoulder_left  = 3;
arm_left       = 4;
forearm_left   = 5;
hand_left      = 6;
shoulder_right = 7;
arm_right      = 8;
forearm_right  = 9;
hand_right     = 10;
abpelv         = 11;
thigh_left     = 12;
leg_left       = 13;
foot_left      = 14;
thigh_right    = 15;
leg_right      = 16;
foot_right     = 17;

% SEGMENTS
person.segment(1).name = 'abdominal-thoracic';
person.segment(2).name = 'head-neck';
person.segment(3).name = 'left shoulder';
person.segment(4).name = 'left arm';
person.segment(5).name = 'left forearm';
person.segment(6).name = 'left hand';
person.segment(7).name = 'right shoulder';
person.segment(8).name = 'right arm';
person.segment(9).name = 'right forearm';
person.segment(10).name = 'right hand';
person.segment(11).name = 'abdominal-pelvic';
person.segment(12).name = 'left thigh';
person.segment(13).name = 'left leg';
person.segment(14).name = 'left foot';
person.segment(15).name = 'right thigh';
person.segment(16).name = 'right leg';
person.segment(17).name = 'right foot';

male = 1;
female = 0;
person.N = 17;

person.solve_ellipse = @(a,u) sqrt(((u/pi).^2)/2-a.^2);
person.resample = @(person,S,x) meas_resample(x,person.meas{S}.length,10,person.segment(S).Nmeas,person.segment(S).Ncalc);
person.plot_points = @(p,varargin) plot3( p(1,:), p(2,:), p(3,:) , varargin{:} );

person.plot = false;
for ii = 1:person.N
  person.segment(ii).plot = false;
end

person.origin{1} = [0;0;0];
for ii = 1:person.N
  person.color{ii} = hsv2rgb( [(ii-1)/person.N , 0.5 , 0.8] );
  person.opacity{ii} = [0.2 0.1]; % face, edge
  person.offset{ii} = [0;0;0];
end