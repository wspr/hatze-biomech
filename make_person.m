%% Make a person

clear all
clc
hfig = figure(1); clf; hold on
set(hfig,'color','white')

male = 1;
female = 0;
i_m = female;

O1 = [0;0;0]; % origin

hatze1 = [...
  101 273 275 264 251 245 262 ... 1:7 ML
  178 ... 8   ???
  084 095 117 149 168 192 190 184 174 174 ... 9:18 AP
  170 ... 19  ???
  452 ... 20  l
  213 ... 21  ???
  ]/1000;
hatze2 = [139 184 214 055]/1000;
hatze11(21) = 0.053;

[calcs,O2] = abdomino_thoracic(O1,i_m,hatze1,hatze2,hatze11);

head_width  = 0.139;
head_depth  = 0.184;
head_height = 0.214;
neck_height = 0.055;

[calcs] = head_neck(O2,i_m,hatze1,head_width,head_depth,head_height,neck_height);

O4 = O2+[0.3;0;0]; % origin of left arm

%left_arm
left_arm_diameters = [085 085 080 077 073 072 072 072 072 083]/1000; 
left_arm_perimeters = [303 292 280 267 260 256 250 243 235 237]/1000; 
left_arm_length = 294/1000;

[calcs, O5] = arm(O4,i_m,'l',left_arm_diameters,left_arm_perimeters,left_arm_length);

O8 = O2+[-0.3;0;0]; % origin of right arm

%right_arm
right_arm_diameters = [093 091 085 087 085 078 074 070 070 076]/1000; 
right_arm_perimeters = [290 279 268 260 257 251 246 236 223 217]/1000; 
right_arm_length = 291/1000;

[calcs, O9] = arm(O8,i_m,'r',right_arm_diameters,right_arm_perimeters,right_arm_length);

%left_forearm
left_forearm_diameters = [085 084 083 082 075 070 065 062 058 055]/1000; 
left_forearm_perimeters = [235 239 239 233 219 201 183 171 162 160]/1000; 
left_forearm_length = 252/1000;

[calcs, O6] = forearm(O5,i_m,'l',left_forearm_diameters,left_forearm_perimeters,left_forearm_length);

%right_forearm
right_forearm_diameters = [079 083 078 076 069 066 066 061 061 056]/1000; 
right_forearm_perimeters = [231 229 225 221 211 195 184 174 164 161]/1000; 
right_forearm_length = 257/1000;

[calcs, O10] = forearm(O9,i_m,'r',right_forearm_diameters,right_forearm_perimeters,right_forearm_length);

%% pelvis

% pelvis height
pelvis_height = 228/1000;

%left_thigh
left_thigh_diameters = [191 191 190 180 162 146 132 126 119 119]/1000; 
left_thigh_perimeters = [629 604 590 565 536 500 453 418 386 373]/1000; 
left_thigh_length = 439/1000;
left_thigh_length2 = 364/1000;

%right_thigh
right_thigh_diameters = [186 187 189 182 171 176 143 136 130 126]/1000; 
right_thigh_perimeters = [620 602 592 563 528 495 461 417 388 377]/1000; 
right_thigh_length = 445/1000;
right_thigh_length2 = 353/1000;

h_l = 0.3*pelvis_height;
h_r = 0.3*pelvis_height; % not sure how these could ever be different?

[calcs, O12, O15] = abdomino_pelvic(O1,i_m,h_l,h_r,...
  left_thigh_diameters, left_thigh_perimeters,...
  right_thigh_diameters,right_thigh_perimeters);

% thighs
[calcs, O13] = thigh(O12,i_m,'l',left_thigh_diameters,left_thigh_perimeters,left_thigh_length,h_l);
[calcs, O16] = thigh(O15,i_m,'r',right_thigh_diameters,right_thigh_perimeters,right_thigh_length,h_r);

%%

%left_leg
left_leg_diameters = [110 103 109 116 117 105 091 079 065 064]/1000; 
left_leg_perimeters = [353 335 362 375 370 338 293 264 238 238]/1000; 
left_leg_length = 406/1000;
left_ankle_size = 0.032;

[calcs, O14] = leg(O13,i_m,'l',left_leg_diameters,left_leg_perimeters,left_leg_length,left_ankle_size);


%right_leg
right_leg_diameters = [109 103 109 116 117 106 095 079 066 062]/1000; 
right_leg_perimeters = [351 332 353 375 374 351 309 279 245 239]/1000; 
right_leg_length = 419/1000;
right_ankle_size = 0.034;

[calcs, O17] = leg(O16,i_m,'r',right_leg_diameters,right_leg_perimeters,right_leg_length,right_ankle_size);

%left_foot
left_foot_ankle_length=[053]/1000;
left_foot_toes_length=[158]/1000;
left_foot_heel_length=[072]/1000;
left_foot_upper_height=[099]/1000;
left_foot_lower_height= [035]/1000;
left_foot_length=[216]/1000;

[calcs]= foot(O14,i_m,'l',left_foot_ankle_length,left_foot_toes_length,left_foot_heel_length,left_foot_upper_height,left_foot_lower_height,left_foot_length);

%right_foot
right_foot_ankle_length=[058]/1000;
right_foot_toes_length=[147]/1000;
right_foot_heel_length=[077]/1000;
right_foot_upper_height=[098]/1000;
right_foot_lower_height= [035]/1000;
right_foot_length=[212]/1000;

[calcs]= foot(O17,i_m,'r',right_foot_ankle_length,right_foot_toes_length,right_foot_heel_length,right_foot_upper_height,right_foot_lower_height,right_foot_length);


if false
plot_coord(O1,'index','1''');
plot_coord(O2,'index','2''');
plot_coord(O4,'index','4''');
plot_coord(O5,'index','5''');
plot_coord(O6,'index','6''');
plot_coord(O12,'index','12''');
plot_coord(O13,'index','13''');
plot_coord(O14,'index','14''');
plot_coord(O16,'index','16''');
plot_coord(O17,'index','17''');
end

axis equal
view(153,23)
axis off
zoom(0.4)