%% PERSON_GENERATE
%
% Main function to calculate and plot Hatze's body segment parameter model

%% Inputs and options
function result = person_generate(person,varargin)

ip = inputParser;
ip.StructExpand = false;
addOptional(ip,'person',struct());
addParamValue(ip,'data','');
addParamValue(ip,'plot','');
addParamValue(ip,'precision','');
parse(ip,person,varargin{:})

%% Main body
%
% * Initialise, if necessary (no person structure input)
% * Load data, if specified (otherwise assume it has already happened)
% * Translate input pose into a form we can process.
% * Then run each calculation function in a convenient loop.

if strcmp(ip.Results.precision,'hatze')
  person.precise = false;
else
  person.precise = true;
end

if isempty(fieldnames(ip.Results.person))
  person = person_initialise();
end

if ~isempty(ip.Results.data)
  person = person_data(person,ip.Results.data);
end

person = person_copy_angles(person);

if ~isempty(ip.Results.plot)
  person.plot = ip.Results.plot;
end

for ii = 1:person.N
  person = person.segment(ii).setup_fn(person,ii);
end

for ii = 1:person.N
  person.segment(ii).Gcentroid = ...
    person.segment(ii).origin + person.segment(ii).Rglobal*person.segment(ii).centroid;
end

if nargout > 0
  result = person;
end

end





%% Subfunction: Initialisation and helper variables

function person = person_initialise()

person.N = 17;

%% Segment names

person.segment( 1).name = 'Abdominal-thoracic';
person.segment( 2).name = 'Head-neck';
person.segment( 3).name = 'Left shoulder';
person.segment( 4).name = 'Left arm';
person.segment( 5).name = 'Left forearm';
person.segment( 6).name = 'Left hand';
person.segment( 7).name = 'Right shoulder';
person.segment( 8).name = 'Right arm';
person.segment( 9).name = 'Right forearm';
person.segment(10).name = 'Right hand';
person.segment(11).name = 'Abdominal-pelvic';
person.segment(12).name = 'Left thigh';
person.segment(13).name = 'Left leg';
person.segment(14).name = 'Left foot';
person.segment(15).name = 'Right thigh';
person.segment(16).name = 'Right leg';
person.segment(17).name = 'Right foot';

%% Setup functions

person.segment( 1).setup_fn = @segment_abdomino_thoracic;
person.segment( 2).setup_fn = @segment_head_neck;
person.segment( 3).setup_fn = @segment_shoulder;
person.segment( 4).setup_fn = @segment_arm;
person.segment( 5).setup_fn = @segment_forearm;
person.segment( 6).setup_fn = @segment_hand;
person.segment( 7).setup_fn = @segment_shoulder;
person.segment( 8).setup_fn = @segment_arm;
person.segment( 9).setup_fn = @segment_forearm;
person.segment(10).setup_fn = @segment_hand;
person.segment(11).setup_fn = @segment_abdomino_pelvic;
person.segment(12).setup_fn = @segment_thigh;
person.segment(13).setup_fn = @segment_leg;
person.segment(14).setup_fn = @segment_foot;
person.segment(15).setup_fn = @segment_thigh;
person.segment(16).setup_fn = @segment_leg;
person.segment(17).setup_fn = @segment_foot;

%% Tree

person.segment( 1).prior = [];
person.segment( 2).prior =  1;
person.segment( 3).prior =  1;
person.segment( 4).prior =  3;
person.segment( 5).prior =  4;
person.segment( 6).prior =  5;
person.segment( 7).prior =  1;
person.segment( 8).prior =  7;
person.segment( 9).prior =  8;
person.segment(10).prior =  9;
person.segment(11).prior =  1;
person.segment(12).prior = 11;
person.segment(13).prior = 12;
person.segment(14).prior = 13;
person.segment(15).prior = 11;
person.segment(16).prior = 15;
person.segment(17).prior = 16;

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
  90; 0;      ... 35, 36     : left foot
  0; 0; 0;   ... 37, 38, 39 : right thigh
  0;         ... 40         : right knee
  90; 0;      ... 41, 42     : right foot
];

%% Inline functions

person.solve_ellipse = @solve_ellipse;
person.resample = @(person,S,x) measurement_resample(x,person.meas{S}.length,10,person.segment(S).Nmeas,person.segment(S).Ncalc);
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

for ii = 1:person.N
  person.segment(ii).colour = hsv2rgb( [(ii-1)/person.N , 0.5 , 0.8] );
  person.segment(ii).opacity = [0.2 0.1]; % face, edge
  person.segment(ii).offset = [0;0;0];
end

end


%% Solve ellipse 
function b = solve_ellipse(a,u)
%gives the second semi-axis of an ellipse, b, given the first semi-axis, a, and
%the perimeter, u. From Biomlib Release 8101 (fortran corrections)

%d=(0.25*u./a);
%ind_small = d<1;
%ind_mid = 1 <= d < pi/2;
%ind_large = d >= pi/2;

%if any(ind_small)
%    aa=a(ind_small);
%    b(ind_small)=0.0001.*aa;
%end

%if any(ind_mid)
%  aa = a(ind_mid);
%  uu = u(ind_mid);
%  b(ind_mid) = aa.*exp((log((uu-4.*aa)/(2*aa*(pi-2))))/1.435);
%end
 
%if any(ind_large)
%  aa = a(ind_large);
%  uu = u(ind_large);
%  b(ind_large)= sqrt(abs((uu./(pi*sqrt(2))).^2-aa.^2));
%end
b=sqrt(abs((u/pi).^2/2-a.^2));
end


%% Subfunction: Read and setup data

function person = person_data(person,filename)

%% Raw measurements

fileID = fopen(filename);
if fileID == -1
  error(['File "',filename,'" not found.'])
end
C = textscan(fileID, '%f','CommentStyle','#');
fclose(fileID);

meas = transpose(C{1}); % rows are easier

person.sex    = double(meas(1));
person.age    = double(meas(2));
person.height = double(meas(3));
person.weight = double(meas(4));
person.scale  = double(meas(5));

meas([1 2 3 4 5]) = [];
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

person.segment( 1).Nmeas = 10;
person.segment( 1).Ncalc = 10;
person.segment( 4).Nmeas = 10;
person.segment( 4).Ncalc = 10;
person.segment( 8).Nmeas = 10;
person.segment( 8).Ncalc = 10;
person.segment( 5).Nmeas = 10;
person.segment( 5).Ncalc = 10;
person.segment( 9).Nmeas = 10;
person.segment( 9).Ncalc = 10;
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
person.meas{1}.depths = person.meas{1}.all(10:19);
person.meas{1}.length = person.meas{1}.all(20);

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

person.meas{11}.diam   = person.meas{11}.all([1:8, 8, 8]);
person.meas{11}.perim  = person.meas{11}.all(11:17);
person.meas{11}.length = person.meas{11}.all(10);

person.meas{12}.diam        = person.meas{12}.all(1:10);
person.meas{12}.perim       = person.meas{12}.all(11:20);
person.meas{12}.length_long = person.meas{12}.all(21);
person.meas{12}.length      = person.meas{12}.all(22);

person.meas{13}.diam   = person.meas{13}.all(1:10);
person.meas{13}.perim  = person.meas{13}.all(11:20);
person.meas{13}.length = person.meas{13}.all(21);
person.meas{13}.ankle  = person.meas{13}.all(22);

person.meas{15}.diam        = person.meas{15}.all(1:10);
person.meas{15}.perim       = person.meas{15}.all(11:20);
person.meas{15}.length_long = person.meas{15}.all(21);
person.meas{15}.length      = person.meas{15}.all(22);

person.meas{16}.diam   = person.meas{16}.all(1:10);
person.meas{16}.perim  = person.meas{16}.all(11:20);
person.meas{16}.length = person.meas{16}.all(21);
person.meas{16}.ankle  = person.meas{16}.all(22);

%% Densities

person.nu = person.meas{11}.all(8)/person.meas{11}.all(9)-1;

person.density.thoracic_wall = @(i_m) 1080+60*i_m;
person.density.abdomen       = @(i_m) 1000+30*i_m;
person.density.lungs         = @(i_m) 300;
person.density.breasts       = @(i_m) 980; 

person.density.head = @(i_m) 1120;
person.density.neck = @(i_m) 1040;

person.density.shoulder_lateral = @(i_m) 1030+20*i_m;
person.density.shoulder_medial  = @(i_m) 1030+20*i_m;
person.density.shoulder_cutout  = @(i_m) 1030+20*i_m;

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

person.density.hand = 1110;


person.density.penis      = 1000;
person.density.thigh_head = @(i_m,nu) 1020 + 20*i_m + 30/((1+2*nu)^2);
person.density.lower_back = @(i_m) 1090 + 30*i_m;
person.density.posterior  = @(i_m) 1020 + 30*i_m;
person.density.stomach    = @(i_m) 1000 + 40*i_m;
person.density.pelvis     = @(i_m) 1020 + 30*i_m;
person.density.buttocks   = @(i_m,nu) 960 + 60*i_m + 30/(1+4*nu)^3;

for ii = 1:3
  person.density.thigh{ii} = @(i_m,nu) 1000+(30+10*(ii-2))/((1+2*nu)^2)+20*i_m;
end
for ii = 4:9
  person.density.thigh{ii} = @(i_m,nu) 1030+10*i_m;
end
  person.density.thigh{10} = @(i_m,nu) 1490+10*i_m;

for ii = 1:3
  person.density.leg{ii} = @(i_m) 1060+90/4*(ii-3)^2;
end
for ii = 4:7
  person.density.leg{ii} = @(i_m) 1060;
end
for ii = 8:10
  person.density.leg{ii} = @(i_m) 1060+130/3*(ii-7);
end

person.density.ankle = 1200;
person.density.foot = @(n) 1480*(1-(0.0001)*(n.^2)*(1-1100/1480));
person.density.heel = 980;
person.density.sole = 1000;

end




%% Subfunction: Copy angles from input to pose
%
% This code translates generalised coordinates into local coordinate system
% rotations for each segment.

function person = person_copy_angles(person)

person.segment( 1).angle = person.q(4:6);
person.segment( 2).angle = person.q(7:9);
person.segment( 3).angle = [person.q(10:11); 0];
person.segment( 4).angle = person.q(12:14);
person.segment( 5).angle = [person.q(15); 0; person.q(16)];
person.segment( 6).angle = [person.q(17:18); 0];
person.segment( 7).angle = [person.q(19:20); 0];
person.segment( 8).angle = person.q(21:23);
person.segment( 9).angle = [person.q(24); 0; person.q(25)];
person.segment(10).angle = [person.q(26:27); 0];
person.segment(11).angle = person.q(28:30);
person.segment(12).angle = person.q(31:33);
person.segment(13).angle = [person.q(34); 0; 0];
person.segment(14).angle = [person.q(35); 0; person.q(36)];
person.segment(15).angle = person.q(37:39);
person.segment(16).angle = [person.q(40); 0; 0];
person.segment(17).angle = [person.q(41); 0; person.q(42)];

for S = 1:person.N
  person.segment(S).Rlocal  = person.cardan_rotation(person.segment(S).angle);
end

person.segment(1).Rglobal = person.segment(1).Rlocal;
for S = 2:person.N
  person.segment(S).Rglobal = ...
    person.segment( person.segment(S).prior ).Rglobal * person.segment(S).Rlocal;
end

end

function person = calculate_centroids(person)

for S = 1:person.N
  person.segment(S).CoM = person.origin{S} + person.segment(S).Rglobal*person.segment(S).centroid;
end

end


function calc = measurement_resample(orig_meas,l,Norig,Nmeas,Ncalc)

orig_pos = linspace(0,l,Norig);
meas_pos = linspace(0,l,Nmeas);
calc_pos = linspace(0,l,Ncalc);

meas = interp1(orig_pos,orig_meas,meas_pos,'cubic');
calc = interp1(meas_pos,meas,calc_pos,'cubic');

end

