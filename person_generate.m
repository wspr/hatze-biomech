

person = segment_abdomino_thoracic(person,1);
person = segment_head_neck(person,2);
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
person = segment_foot(person,14,'l');

person = segment_thigh(person,15,'r',h_r);
person = segment_leg(person,16,'r',right_leg_diameters,right_leg_perimeters,right_leg_length,right_ankle_size);
person = segment_foot(person,17,'r');

