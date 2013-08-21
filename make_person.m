%% Make a person

clear all
clc

male = 1;
female = 0;
i_m = female;

% SEGMENTS
%  1 abdom-thor
%  2 head-neck
%  3 left shoulder
%  4 left arm
%  5 left forearm
%  6 left hand
%  7 right shoulder
%  8 right arm
%  9 right forearm
% 10 right hand
% 11 abdom-pelvic
% 12 left thigh
% 13 left leg
% 14 left foot
% 15 right thigh
% 16 right left
% 17 right foot

person.N = 17;

person.sex = female;
person.origin{1} = [0;0;0];
for ii = 1:person.N
  person.color{ii} = hsv2rgb( [(ii-1)/person.N , 0.5 , 0.8] );
  person.opacity{ii} = [0.2 0.1]; % face, edge
  person.offset{ii} = [0;0;0];
end

%% offsets for visualisation

person.offset{1} =  [0; 0;  0];
person.offset{2} =  [0; 0; 50];
person.offset{3} =  [-100; 0; 50];
person.offset{4} =  [0; 0; -100];
person.offset{5} =  [0; 0; -50];
person.offset{6} =  [0; 0; -50];
person.offset{7} =  [100; 0; 50];
person.offset{8} =  [0; 0; -100];
person.offset{9} =  [0; 0; -50];
person.offset{10} = [0; 0; -50];

person.offset{11} = [0; 0; -50]; % pelvis
person.offset{12} = [0; 0; -100];
person.offset{13} = [0; 0; -50];
person.offset{14} = [0; 0; -50];
person.offset{15} = [0; 0; -100];
person.offset{16} = [0; 0; -50];
person.offset{17} = [0; 0; -50];

%%

person_data

%%

hfig = figure(1); clf; hold on
set(hfig,'color','white')

person = segment_abdomino_thoracic(person,1,hatze1,hatze2,hatze11);
person = segment_head_neck(person,2,hatze1,head_width,head_depth,head_height,neck_height);
person = segment_abdomino_pelvic(person,11,...
  pelvis_widths,pelvis_perimeters,pelvis_meas,...
  h_l,h_r);

person = segment_shoulder(person,3,'l',left_arm_diameters);
person = segment_arm(person,4,'l',left_arm_diameters,left_arm_perimeters,left_arm_length);
person = segment_forearm(person,5,'l',left_forearm_diameters,left_forearm_perimeters,left_forearm_length);
person = segment_hand(person,6);

person = segment_shoulder(person,7,'r',right_arm_diameters);
person = segment_arm(person,8,'r',right_arm_diameters,right_arm_perimeters,right_arm_length);
person = segment_forearm(person,9,'r',right_forearm_diameters,right_forearm_perimeters,right_forearm_length);
person = segment_hand(person,10);

person = segment_thigh(person,12,'l',h_l);
person = segment_leg(person,13,'l',left_leg_diameters,left_leg_perimeters,left_leg_length,left_ankle_size);
person = segment_foot(person,14,'l',left_foot_ankle_length,left_foot_toes_length,left_foot_heel_length,left_foot_upper_height,left_foot_lower_height,left_foot_length);

person = segment_thigh(person,15,'r',h_r);
person = segment_leg(person,16,'r',right_leg_diameters,right_leg_perimeters,right_leg_length,right_ankle_size);
person = segment_foot(person,17,'r',right_foot_ankle_length,right_foot_toes_length,right_foot_heel_length,right_foot_upper_height,right_foot_lower_height,right_foot_length);

axis equal
view(153,23)
axis off
zoom(2)

%%

ind = 1:person.N;
ind([3 7]) = []; % repeated
for ii = ind
%  plot_coord(person.origin{ii},'index',[num2str(ii),'''']);
end

%%

plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',2)
