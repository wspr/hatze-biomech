%% Initialisation and helper variables
%
% This code should be run before any other person_*.m files.

%% Logical indexing
%
% Not sure if I want to use these or not.

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

male = 1;
female = 0;


%% Segment names

person.N = 17;

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

%% Setup functions

person.segment(1).setup_fn = @segment_abdomino_thoracic;
person.segment(2).setup_fn = @segment_head_neck;
person.segment(3).setup_fn = @segment_shoulder;
person.segment(4).setup_fn = @segment_arm;
person.segment(5).setup_fn = @segment_forearm;
person.segment(6).setup_fn = @segment_hand;
person.segment(7).setup_fn = @segment_shoulder;
person.segment(8).setup_fn = @segment_arm;
person.segment(9).setup_fn = @segment_forearm;
person.segment(10).setup_fn = @segment_hand;
person.segment(11).setup_fn = @segment_abdomino_pelvic;
person.segment(12).setup_fn = @segment_thigh;
person.segment(13).setup_fn = @segment_leg;
person.segment(14).setup_fn = @segment_foot;
person.segment(15).setup_fn = @segment_thigh;
person.segment(16).setup_fn = @segment_leg;
person.segment(17).setup_fn = @segment_foot;


%% Joint angles

person.q = [ ...
  0; 0; 0;   ...  1,  2,  3 : global position 
  0; 0; 0;   ...  4,  5,  6 : orientation of segment 1 (abdomen-thorax)
  0; 0; 0;   ...  7,  8,  9 : orientation of segment 2 (head-neck)
  0; -90;    ... 10, 11     : left shoulder
  0;  90; 0; ... 12, 13, 14 : left arm
  0;  90;    ... 15, 16     : left forearm
  0; 0;      ... 17, 18     : left hand
  0;  90;    ... 19, 20     : right shoulder
  0; -90; 0; ... 21, 22, 23 : right arm
  0; -90;    ... 24, 25     : right forearm
  0; 0;      ... 26, 27     : right hand
  0; 0; 0;   ... 28, 29, 30 : abdomen-pelvis
  0; 0; 0;   ... 31, 32, 33 : left thigh
  0;         ... 34         : left knee
  0; 0;      ... 35, 36     : left foot
  0; 0; 0;   ... 37, 38, 39 : right thigh
  0;         ... 40         : right knee
  0; 0;      ... 41, 42     : right foot
];

%% Inline functions

person.solve_ellipse = @(a,u) sqrt(((u/pi).^2)/2-a.^2);
person.resample = @(person,S,x) meas_resample(x,person.meas{S}.length,10,person.segment(S).Nmeas,person.segment(S).Ncalc);
person.plot_points = @(p,varargin) plot3( p(1,:), p(2,:), p(3,:) , varargin{:} );
person.cardan_rotation = @(a) ...
  [1 0 0;0 cosd(a(1)) -sind(a(1));0 sind(a(1)) cosd(a(1))]*...
  [cosd(a(2)) 0 sind(a(2));0 1 0;-sind(a(2)) 0 cosd(a(2))]*...
  [cosd(a(3)) -sind(a(3)) 0;sind(a(3)) cosd(a(3)) 0;0 0 1];

%% Plotting defaults

person.plot = false;
for ii = 1:person.N
  person.segment(ii).plot = false;
end

person.origin{1} = person.q(1:3);
person.segment(1).frame = [person.cardan_rotation(person.q(4:6)) person.q(1:3);0 0 0 1];

for ii = 1:person.N
  person.color{ii} = hsv2rgb( [(ii-1)/person.N , 0.5 , 0.8] );
  person.opacity{ii} = [0.2 0.1]; % face, edge
  person.offset{ii} = [0;0;0];
end

