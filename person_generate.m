

person = segment_abdomino_thoracic(person,1);
person = segment_head_neck(person,2);
person = segment_abdomino_pelvic(person,11,...
  pelvis_widths,pelvis_perimeters,pelvis_meas,...
  h_l,h_r);

% left
person = segment_shoulder(person,3);
person = segment_arm(person,4);
person = segment_forearm(person,5);
person = segment_hand(person,6);

% right
person = segment_shoulder(person,7);
person = segment_arm(person,8);
person = segment_forearm(person,9);
person = segment_hand(person,10);

% left
person = segment_thigh(person,12);
person = segment_leg(person,13);
person = segment_foot(person,14);

% right
person = segment_thigh(person,15);
person = segment_leg(person,16);
person = segment_foot(person,17);

