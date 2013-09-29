%% Driver file to experiment with movement

clear all
clc

person = person_generate('data','hatze_meas.txt');

person.q = [ ...
  0; 0; 0;   ...  1,  2,  3 : global *position*; rest are angles:
  20; 0; 0;   ...  4,  5,  6 : abdomen-thorax
  45; 0; 0;   ...  7,  8,  9 : head-neck
  45; -135;    ... 10, 11     : left shoulder
  0;  -90; 90; ... 12, 13, 14 : left arm
  90;  0;    ... 15, 16     : left forearm
  0; 0;      ... 17, 18     : left hand
  0;  55;    ... 19, 20     : right shoulder
  45; 90; 0; ... 21, 22, 23 : right arm
  45; -90;    ... 24, 25     : right forearm
  0; 0;      ... 26, 27     : right hand
  30; 0; 0;   ... 28, 29, 30 : abdomen-pelvis
  90; 0; 0;   ... 31, 32, 33 : left thigh
  -90;         ... 34         : left knee
  90; 0;      ... 35, 36     : left foot
  60; -45; 0;   ... 37, 38, 39 : right thigh
  -90;         ... 40         : right knee
  90; 0;      ... 41, 42     : right foot
];

if false
person.q(4)  =  person.segment( 1).theta;
person.q(7)  = -person.segment( 1).theta;
person.q(10) = -person.segment( 1).theta;
person.q(18) =  person.segment( 6).theta;
person.q(19) = -person.segment( 1).theta;
person.q(27) =  person.segment(10).theta;
person.q(28) =  person.segment(11).theta - person.segment(1).theta;
person.q(31) = -person.segment(11).theta;
person.q(35) =  person.segment(14).theta + 90;
person.q(37) = -person.segment(11).theta;
person.q(41) =  person.segment(17).theta + 90;
end

figure(111); clf; hold on

person.plot = true; % all segments
%person.segment(4).plot = true;
%person.segment(8).plot = true;
person = person_generate(person);

person.plot_points([person.segment(3:4).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(7:8).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([1,3,2]).origin], 'k.-', 'markersize', 20,'linewidth',2)

axis equal
view([140 20])
pbaspect([1 1 1])

%%

person.plot_points([person.segment(3:6).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(7:10).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(11:14).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([11,15:17]).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([1,3,2]).origin], 'k.-', 'markersize', 20,'linewidth',2)

ind = 1:person.N;
for ii = ind
  plot_coord(person.segment(ii).origin,'index',[num2str(ii),''''],'rotate',person.segment(ii).Rglobal,'length',0.07);
end

for ii = 1:person.N  person.plot_points(person.segment(ii).Gcentroid, 'r.', 'markersize', 30)
end

axis equal
view([140 20])
pbaspect([1 1 1])
