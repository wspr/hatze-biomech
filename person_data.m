function person = person_data(person,filename)

%% Raw measurements

fileID = fopen(filename);
if fileID == -1
  error(['File "',filename,'" not found.'])
end
C = textscan(fileID, '%d','CommentStyle','#');
fclose(fileID);

meas = transpose(C{1}); % rows are easier
omeas = meas;

person.sex = double(meas(1));
person.scale = double(meas(2));

meas([1 2]) = [];
meas = double(meas)/person.scale;

% order of segments in the measurement file:
seg   = [ 1 2 3 7  4  8  5  9 6 10 11 12 15 13 16 14 17];

% number of measurements per segment:
count = [21 4 4 4 21 21 21 21 2  2 21 22 22 22 22  6  6];

for S = 1:17
  
  person.meas{seg(S)}.all = meas(1:count(S));
  meas(1:count(S)) = [];
  
end

%% Index counts
%
% Pretty sure I can do better than this!

person.segment(1).Nmeas = 10;
person.segment(1).Ncalc = 10;
person.segment(4).Nmeas = 10;
person.segment(4).Ncalc = 10;
person.segment(8).Nmeas = 10;
person.segment(8).Ncalc = 10;
person.segment(5).Nmeas = 10;
person.segment(5).Ncalc = 10;
person.segment(9).Nmeas = 10;
person.segment(9).Ncalc = 10;
person.segment(11).Nmeas = 10;
person.segment(11).Ncalc = 10;
person.segment(12).Nmeas = 10;
person.segment(12).Ncalc = 10;
person.segment(15).Nmeas = 10;
person.segment(15).Ncalc = 10;
person.segment(13).Nmeas = 10;
person.segment(13).Ncalc = 10;
person.segment(16).Nmeas = 10;
person.segment(16).Ncalc = 10;

%% Organising measurements

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

person.meas{11}.diam   = person.meas{11}.all(1:10);
person.meas{11}.perim  = person.meas{11}.all(11:17);
person.meas{11}.length = person.meas{11}.all(20);

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

%% Densities

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
for n = 2:9
  person.density.arm{n} = @(i_m) 1058+20*i_m;
end
person.density.arm{10} = @(i_m) 1080+20*i_m;
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

end
