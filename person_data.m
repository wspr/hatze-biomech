
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
person.segment(1).N = 10;

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
person.segment(4).N = 10;

person.meas{8}.all = [...
    093 091 085 087 085 078 074 070 070 076 ... diameters
    290 279 268 260 257 251 246 236 223 217 ... perimeters
    291 ... length
  ];
person.segment(8).N = 10;

% 5 & 9: forearms

person.meas{5}.all = [...
    085 084 083 082 075 070 065 062 058 055 ... diameters
    235 239 239 233 219 201 183 171 162 160 ... perimeters
    252 ... length
  ];
person.segment(5).N = 10;

person.meas{9}.all = [...
    079 083 078 076 069 066 066 061 061 056 ... diameters
    231 229 225 221 211 195 184 174 164 161 ... perimeters
    257 ... length
  ];
person.segment(9).N = 10;

% 6 & 10: hands
person.meas{6}.all  = [43 77];
person.meas{10}.all = [45 76];

% 11: pelvis

person.meas{11}.all = [...
    287 301 318 331 341 344 359 382 326 240 ...
    767 797 844 887 940 954 975 ...
    097 081 228 053 ...
  ];
person.segment(11).N = 10;

% 12 & 15 thighs

person.meas{12}.all = [...
    191 191 190 180 162 146 132 126 119 119 ... diam
    629 604 590 565 536 500 453 418 386 373 ... perim
    439 ... ???
    364 ... length
];
person.segment(12).N = 10;

person.meas{15}.all = [...
    186 187 189 182 171 176 143 136 130 126 ... diam
    620 602 592 563 528 495 461 417 388 377 ... perim
    445 ... ???
    353 ... length
];
person.segment(15).N = 10;

% 13 & 16: legs
person.meas{13}.all = [...
    110 103 109 116 117 105 091 079 065 064 ... diam
    353 335 362 375 370 338 293 264 238 238 ... perim
    406 ... length
    032 ... ankle size
];
person.segment(13).N = 10;

person.meas{16}.all = [...
    109 103 109 116 117 106 095 079 066 062 ... diam
    351 332 353 375 374 351 309 279 245 239 ... perim
    419 ... length
    034 ... ankle size
];
person.segment(16).N = 10;

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

for ii = 1:person.N
  if strcmp(person.units,'mm')
    person.meas{ii}.all = person.meas{ii}.all/1000;
  end
end

person.meas{1}.widths = [person.meas{1}.all(1) NaN NaN NaN person.meas{1}.all(2:7)];
person.meas{1}.depths = person.meas{1}.all(9:18);

person.meas{4}.diam   = person.meas{4}.all(1:10); 
person.meas{4}.perim  = person.meas{4}.all(11:20); 
person.meas{4}.length = person.meas{4}.all(21);

person.meas{8}.diam   = person.meas{8}.all(1:10); 
person.meas{8}.perim  = person.meas{8}.all(11:20); 
person.meas{8}.length = person.meas{8}.all(21);

person.meas{5}.diam   = person.meas{5}.all(1:10); 
person.meas{5}.perim  = person.meas{5}.all(11:20); 
person.meas{5}.length = person.meas{5}.all(21);

person.meas{9}.diam   = person.meas{9}.all(1:10); 
person.meas{9}.perim  = person.meas{9}.all(11:20); 
person.meas{9}.length = person.meas{9}.all(21);

pelvis_widths     = person.meas{11}.all(1:10);
pelvis_perimeters = person.meas{11}.all(11:17);
pelvis_meas       = person.meas{11}.all(18:21);
pelvis_height     = person.meas{11}.all(20);

h_l = 0.3*pelvis_height;
h_r = 0.3*pelvis_height; % not sure how these could ever be different?

person.meas{12}.diam   = person.meas{12}.all(1:10);
person.meas{12}.perim  = person.meas{12}.all(11:20);
person.meas{12}.length_long = person.meas{12}.all(21);
person.meas{12}.length = person.meas{12}.all(22);

person.meas{13}.diam   = person.meas{13}.all(1:10);
person.meas{13}.perim  = person.meas{13}.all(11:20);
person.meas{13}.length = person.meas{13}.all(21);
person.meas{13}.ankle  = person.meas{13}.all(22);

person.meas{15}.diam   = person.meas{15}.all(1:10);
person.meas{15}.perim  = person.meas{15}.all(11:20);
person.meas{15}.length_long = person.meas{15}.all(21);
person.meas{15}.length = person.meas{15}.all(22);

person.meas{16}.diam   = person.meas{16}.all(1:10);
person.meas{16}.perim  = person.meas{16}.all(11:20);
person.meas{16}.length = person.meas{16}.all(21);
person.meas{16}.ankle  = person.meas{16}.all(22);

%%

person.nu = 0.1;

person.density.thoracic_wall = @(i_m) 1080+60*i_m;
person.density.abdomen       = @(i_m) 1000+30*i_m;
person.density.lungs         = @(i_m) 300;
person.density.breasts       = @(i_m) 1200; % (? see A2.94)

person.density.head = @(i_m) 1120;
person.density.neck = @(i_m) 1040;

person.density.shoulder_lateral = @(i_m) 1030+20*i_m;
person.density.shoulder_medial  = @(i_m) 1030+20*i_m;
person.density.shoulder_cutout  = @(i_m) 1030+20*i_m;

person.density.penis      = 1000;
person.density.thigh_head = @(i_m,nu) 1020 + 20*i_m + 30/((1+2*nu)^2); 
person.density.lower_back = @(i_m) 1090 + 30*i_m;
person.density.posterior  = @(i_m) 1020 + 30*i_m;
person.density.stomach    = @(i_m) 1000 + 40*i_m;
person.density.pelvis     = @(i_m) 1020 + 30*i_m;
person.density.buttocks   = @(i_m,nu) 960 + 60*i_m + 30/(1+4*nu)^3;

person.density.arm{1} = @(i_m) 1060+40*i_m;
for n = 2:(person.segment(4).N-1)
  person.density.arm{n} = @(i_m) 1058+20*i_m;
end
person.density.arm{person.segment(4).N} = @(i_m) 1080+20*i_m;
person.density.humerous = @(i_m) 1080+20*i_m;

for n = 1:2
  person.density.forearm{n} = @(i_m) (1160-60*n)*(1+0.0213*i_m);
end
for n = 3:8
  person.density.forearm{n} = @(i_m) (1034.2+2.86*n)*(1+0.0213*i_m);
end
for n = 9:10
  person.density.forearm{n} = @(i_m) 1204.29*(1+0.0213*i_m);
end

for ii = 1:3
  person.density.thigh{ii} = @(i_m,nu) 1000+(30+10*(ii-2))/((1+2*nu)^2)+20*i_m;
end
for ii = 4:9
  person.density.thigh{ii} = @(i_m,nu) 1030+10*i_m;
end
person.density.thigh{10} = @(i_m,nu) 1490+10*i_m;

for ii = 1:3
  person.density.leg{ii} = @(i_m) 1060+22.5*(ii-3)^2;
end
for ii = 4:7
  person.density.leg{ii} = @(i_m) 1060;
end
for ii = 8:10
  person.density.leg{ii} = @(i_m) 1060+43.33*(ii-7);
end
person.density.ankle = 1200;
person.density.foot = @(n) 1480*(1-(0.0001)*(n.^2)*(1-1100/1480));


%%

% left arm
person.segment(4).mass_hatze = 1.616;
person.segment(4).volume_hatze = 1.713;
person.segment(4).centroid_hatze = [0; 0; -0.131];
person.segment(4).Minertia_hatze = [0.013677 0.013435 0.001528];

% left forearm
person.segment(5).mass_hatze = 0.835;
person.segment(5).volume_hatze = 0.896;
person.segment(5).centroid_hatze = [0; 0; -0.108];
person.segment(5).Minertia_hatze = [0.004596 0.004712 0.000523];

% right arm
person.segment(8).mass_hatze = 1.505;
person.segment(8).volume_hatze = 1.595;
person.segment(8).centroid_hatze = [0; 0; -0.129];
person.segment(8).Minertia_hatze = [0.011917 0.011937 0.001325];

% right forearm
person.segment(9).mass_hatze = 0.809;
person.segment(9).volume_hatze = 0.869;
person.segment(9).centroid_hatze = [0; 0; -0.112];
person.segment(9).Minertia_hatze = [0.004765 0.004855 0.000472];

% left thigh
person.segment(12).mass_hatze = 9.166;
person.segment(12).volume_hatze = 9.593;
person.segment(12).centroid_hatze = [0; 0; -0.194];
person.segment(12).Minertia_hatze = [0.156305 0.152523 0.035853];

% left leg
person.segment(13).mass_hatze = 3.310;
person.segment(13).volume_hatze = 3.602;
person.segment(13).centroid_hatze = [0; 0; -0.174];
person.segment(13).Minertia_hatze = [0.045315 0.044973 0.005106];

% right thigh
person.segment(15).mass_hatze = 8.955;
person.segment(15).volume_hatze = 9.375;
person.segment(15).centroid_hatze = [0; 0; -0.194];
person.segment(15).Minertia_hatze = [0.149153 0.146883 0.034605];

% right leg
person.segment(16).mass_hatze = 3.487;
person.segment(16).volume_hatze = 3.794;
person.segment(16).centroid_hatze = [0; 0; -0.183];
person.segment(16).Minertia_hatze = [0.050537 0.050114 0.005438];
