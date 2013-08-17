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

person.meas(1).widths = [101 NaN NaN NaN 273 275 264 251 245 262];
person.meas(1).depths = [084 095 117 149 168 192 190 184 174 174];
hatze1 = [...
  101 273 275 264 251 245 262 ... 1:7 ML
  178 ... 8   ???
  084 095 117 149 168 192 190 184 174 174 ... 9:18 AP
  170 ... 19  ???
  452 ... 20  l
  213 ... 21  ???
  ];
hatze2 = [139 184 214 055];
hatze11(21) = 53;

head_width  = 139;
head_depth  = 184;
head_height = 214;
neck_height = 055;

left_arm_diameters = [085 085 080 077 073 072 072 072 072 083]; 
left_arm_perimeters = [303 292 280 267 260 256 250 243 235 237]; 
left_arm_length = 294;

%left_forearm
left_forearm_diameters = [085 084 083 082 075 070 065 062 058 055]; 
left_forearm_perimeters = [235 239 239 233 219 201 183 171 162 160]; 
left_forearm_length = 252;

%right_arm
right_arm_diameters = [093 091 085 087 085 078 074 070 070 076]; 
right_arm_perimeters = [290 279 268 260 257 251 246 236 223 217]; 
right_arm_length = 291;

%right_forearm
right_forearm_diameters = [079 083 078 076 069 066 066 061 061 056]; 
right_forearm_perimeters = [231 229 225 221 211 195 184 174 164 161]; 
right_forearm_length = 257;

% pelvis
pelvis_widths = [287 301 318 331 341 344 359 382 326 240];
pelvis_perimeters = [767 797 844 887 940 954 975];
pelvis_meas = [097 081 228 053];
pelvis_height = 228;
h_l = 0.3*pelvis_height;
h_r = 0.3*pelvis_height; % not sure how these could ever be different?

%left_thigh
left_thigh_diameters = [191 191 190 180 162 146 132 126 119 119]; 
left_thigh_perimeters = [629 604 590 565 536 500 453 418 386 373]; 
left_thigh_length = 439;
left_thigh_length2 = 364;

%left_leg
left_leg_diameters = [110 103 109 116 117 105 091 079 065 064]; 
left_leg_perimeters = [353 335 362 375 370 338 293 264 238 238]; 
left_leg_length = 406;
left_ankle_size = 032;

%left_foot
left_foot_ankle_length=[053];
left_foot_toes_length=[158];
left_foot_heel_length=[072];
left_foot_upper_height=[099];
left_foot_lower_height= [035];
left_foot_length=[216];

%right_thigh
right_thigh_diameters = [186 187 189 182 171 176 143 136 130 126]; 
right_thigh_perimeters = [620 602 592 563 528 495 461 417 388 377]; 
right_thigh_length = 445;
right_thigh_length2 = 353;

%right_leg
right_leg_diameters = [109 103 109 116 117 106 095 079 066 062]; 
right_leg_perimeters = [351 332 353 375 374 351 309 279 245 239]; 
right_leg_length = 419;
right_ankle_size = 034;

%right_foot
right_foot_ankle_length=[058];
right_foot_toes_length=[147];
right_foot_heel_length=[077];
right_foot_upper_height=[098];
right_foot_lower_height= [035];
right_foot_length=[212];

%%

hfig = figure(1); clf; hold on
set(hfig,'color','white')

[thorax_length,person] = abdomino_thoracic(person,hatze1,hatze2,hatze11);
person = head_neck(person,hatze1,head_width,head_depth,head_height,neck_height);
person = shoulder(person,'l',thorax_length,left_arm_diameters);
person = arm(person,'l',left_arm_diameters,left_arm_perimeters,left_arm_length);
person = forearm(person,'l',left_forearm_diameters,left_forearm_perimeters,left_forearm_length);
person = shoulder(person,'r',thorax_length,right_arm_diameters);
person = arm(person,'r',right_arm_diameters,right_arm_perimeters,right_arm_length);
person = forearm(person,'r',right_forearm_diameters,right_forearm_perimeters,right_forearm_length);

person = abdomino_pelvic(person,...
  pelvis_widths,pelvis_perimeters,pelvis_meas,...
  h_l,h_r,...
  left_thigh_diameters, left_thigh_perimeters,...
  right_thigh_diameters,right_thigh_perimeters);

person = thigh(person,'l',left_thigh_diameters,left_thigh_perimeters,left_thigh_length,h_l);
person = leg(person,'l',left_leg_diameters,left_leg_perimeters,left_leg_length,left_ankle_size);
person = foot(person,'l',left_foot_ankle_length,left_foot_toes_length,left_foot_heel_length,left_foot_upper_height,left_foot_lower_height,left_foot_length);
person = thigh(person,'r',right_thigh_diameters,right_thigh_perimeters,right_thigh_length,h_r);
person = leg(person,'r',right_leg_diameters,right_leg_perimeters,right_leg_length,right_ankle_size);
person = foot(person,'r',right_foot_ankle_length,right_foot_toes_length,right_foot_heel_length,right_foot_upper_height,right_foot_lower_height,right_foot_length);

axis equal
view(153,23)
axis off
zoom(2)

%%

ind = 1:person.N;
ind([3 7]) = []; % repeated
for ii = ind
  plot_coord(person.origin{ii},'index',[num2str(ii),'''']);
end

%%

plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',4)
plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',4)
plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',4)
plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',4)
plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',4)
