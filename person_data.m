
%% Raw measurements

% abdominal-thoracic
person.meas{1}.all = [...
  101 273 275 264 251 245 262 ... 1:7 ML widths
  178 ... 8   ???
  084 095 117 149 168 192 190 184 174 174 ... 9:18 AP depths
  170 ... 19  ???
  452 ... 20  length
  213 ... 21  ???
  ];

% head-neck
person.meas{2}.all = [139 184 214 055]; % width, depth, height, neck height

% 3 & 7: shoulders
person.meas{3}.all = [140 176 090 025]; % d, 2*b, ?, z_h
person.meas{7}.all = [145 182 086 025];

% 4 & 8: arms
person.meas{4}.all = [...
    085 085 080 077 073 072 072 072 072 083 ... diameters
    303 292 280 267 260 256 250 243 235 237 ... perimeters
    294 ... length
  ];

person.meas{8}.all = [...
    093 091 085 087 085 078 074 070 070 076 ... diameters
    290 279 268 260 257 251 246 236 223 217 ... perimeters
    291 ... length
  ];

% 5 & 9: forearms

person.meas{5}.all = [...
    085 084 083 082 075 070 065 062 058 055 ... diameters
    235 239 239 233 219 201 183 171 162 160 ... perimeters
    252 ... length
  ];

person.meas{9}.all = [...
    079 083 078 076 069 066 066 061 061 056 ... diameters
    231 229 225 221 211 195 184 174 164 161 ... perimeters
    257 ... length
  ];

% 6 & 10: hands
person.meas{6}.all  = [43 77];
person.meas{10}.all = [45 76];

% 11: pelvis

person.meas{11}.all = [...
    287 301 318 331 341 344 359 382 326 240 ...
    767 797 844 887 940 954 975 ...
    097 081 228 053 ...
  ];

% 12 & 15 thighs

person.meas{12}.all = [...
    191 191 190 180 162 146 132 126 119 119 ... diam
    629 604 590 565 536 500 453 418 386 373 ... perim
    364 ... length
    439 ... ???
];

person.meas{15}.all = [...
    186 187 189 182 176 171 143 136 130 126 ... diam NB 5/6 swapped from Hatze!
    620 602 592 563 528 495 461 417 388 377 ... perim
    353 ... length
    445 ... ???
];

% 13 & 16: legs
person.meas{13}.all = [...
    110 103 109 116 117 105 091 079 065 064 ... diam
    353 335 362 375 370 338 293 264 238 238 ... perim
    406 ... length
    032 ... ankle size
];

person.meas{16}.all = [...
    109 103 109 116 117 106 095 079 066 062 ... diam
    351 332 353 375 374 351 309 279 245 239 ... perim
    419 ... length
    034 ... ankle size
];

% 14 & 17: feet

person.meas{14}.all = [...
    053 ... ankle length
    158 ... toes length
    072 ... heel length
    035 ... lower height
    099 ... upper+lower height
    216 ... length
];

person.meas{17}.all = [...
    058 ... ankle length
    147 ... toes length
    077 ... heel length
    035 ... lower height
    098 ... upper+lower height
    212 ... length
];

%%

person.meas{1}.widths = [person.meas{1}.all(1) NaN NaN NaN person.meas{1}.all(2:7)];
person.meas{1}.depths = person.meas{1}.all(9:18);

left_arm_diameters  = person.meas{4}.all(1:10); 
left_arm_perimeters = person.meas{4}.all(11:20); 
left_arm_length     = person.meas{4}.all(21);

right_arm_diameters  = person.meas{8}.all(1:10); 
right_arm_perimeters = person.meas{8}.all(11:20); 
right_arm_length     = person.meas{8}.all(21);

left_forearm_diameters  = person.meas{5}.all(1:10); 
left_forearm_perimeters = person.meas{5}.all(11:20); 
left_forearm_length     = person.meas{5}.all(21);

right_forearm_diameters  = person.meas{9}.all(1:10); 
right_forearm_perimeters = person.meas{9}.all(11:20); 
right_forearm_length     = person.meas{9}.all(21);

pelvis_widths     = person.meas{11}.all(1:10);
pelvis_perimeters = person.meas{11}.all(11:17);
pelvis_meas       = person.meas{11}.all(18:21);
pelvis_height     = person.meas{11}.all(20);

h_l = 0.3*pelvis_height;
h_r = 0.3*pelvis_height; % not sure how these could ever be different?

person.meas{12}.diam   = person.meas{12}.all(1:10);
person.meas{12}.perim  = person.meas{12}.all(11:20);
person.meas{12}.length = person.meas{12}.all(21);

person.meas{15}.diam   = person.meas{15}.all(1:10);
person.meas{15}.perim  = person.meas{15}.all(11:20);
person.meas{15}.length = person.meas{15}.all(21);

left_leg_diameters  = person.meas{13}.all(1:10); 
left_leg_perimeters = person.meas{13}.all(11:20); 
left_leg_length     = person.meas{13}.all(21);
left_ankle_size     = person.meas{13}.all(22);

right_leg_diameters  = person.meas{16}.all(1:10); 
right_leg_perimeters = person.meas{16}.all(11:20); 
right_leg_length     = person.meas{16}.all(21);
right_ankle_size     = person.meas{16}.all(22);

