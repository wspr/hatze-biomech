%% Make a person

clear all
clc

% SEGMENTS
person.segment(1).name = 'abdominal-thoracic';
person.segment(2).name = 'head-neck';
person.segment(3).name = 'left shoulder';
person.segment(4).name = 'left arm';
person.segment(5).name = 'left forearm';
person.segment(6).name = 'left hand';
person.segment(7).name = 'right shoulder';
person.segment(8).name = 'right arm';
person.segment(9).name = 'right forearm';
person.segment(10).name = 'right hand';
person.segment(11).name = 'abdominal-pelvic';
person.segment(12).name = 'left thigh';
person.segment(13).name = 'left leg';
person.segment(14).name = 'left foot';
person.segment(15).name = 'right thigh';
person.segment(16).name = 'right left';
person.segment(17).name = 'right foot';

male = 1;
female = 0;
i_m = female;

person.N = 17;

person.sex = female;
person.origin{1} = [0;0;0];
for ii = 1:person.N
  person.color{ii} = hsv2rgb( [(ii-1)/person.N , 0.5 , 0.8] );
  person.opacity{ii} = [0.2 0.1]; % face, edge
  person.offset{ii} = [0;0;0];
end

person_data

count = 0;
for ii = 1:person.N
%    fprintf('Segment: %s\nMeasurements: %i\n',person.segment(ii).name,numel(person.meas{ii}.all))
    fprintf('%s & %i\\\\\n',person.segment(ii).name,numel(person.meas{ii}.all))
    count = count + numel(person.meas{ii}.all);
end
count

%%

figure(1); clf; hold on
set(gcf,'color','white')

person_generate

axis equal
view(153,23)
axis off
zoom(2)

plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',2)
plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',2)

pbaspect([2 1 4])
figuresize(8,12,'centimeters')
matlabfrag('fig/hatze-iso','renderer','opengl')

%%

ind = 1:person.N;
ind([3 7]) = []; % repeated
for ii = ind
%  plot_coord(person.origin{ii},'index',[num2str(ii),'''']);
end

%% FRONT and SIDE

% offsets for visualisations below

person.offset{1} =  [0; 0;  0];
person.offset{2} =  [0; 0; 50];
person.offset{3} =  [-100; 250; 100]; % left shoulder
person.offset{4} =  [-50; 0; -140];
person.offset{5} =  [0; 0; -50];
person.offset{6} =  [0; 60; -50];
person.offset{7} =  [100; -250; 100]; % right shoulder
person.offset{8} =  [50; 0; -140];
person.offset{9} =  [0; 0; -50];
person.offset{10} = [0; -60; -50];

person.offset{11} = [0; 0; -50]; % pelvis
person.offset{12} = [0; 120; -140];
person.offset{13} = [0; 0; -50];
person.offset{14} = [0; 0; -50];
person.offset{15} = [0; -120; -140];
person.offset{16} = [0; 0; -50];
person.offset{17} = [0; 0; -50];

figure(2); clf; hold on
set(gcf,'color','white')

person_generate

axis equal
view(153,23)
axis off
zoom(4)

view([180 0])

pbaspect([1 4 4])
figuresize(8,16,'centimeters')
matlabfrag('fig/hatze-front','renderer','opengl')

figure(3); clf; hold on
set(gcf,'color','white')

person_generate

axis equal
view(153,23)
axis off
zoom(4)

view([90 0])

pbaspect([1 4 4])
figuresize(8,16,'centimeters')
matlabfrag('fig/hatze-side','renderer','opengl')



%% TOP

% offsets for visualisations below

person.offset{1} =  [0; 0;  0];
person.offset{2} =  [0; -250;  0];
person.offset{3} =  [-150; 0;  0]; % left shoulder
person.offset{4} =  [-80; 0;  0];
person.offset{5} =  [-120; 0;  0];
person.offset{6} =  [-120; 0;  0];
person.offset{7} =  [150; 0;  0]; % right shoulder
person.offset{8} =  [80; 0;  0];
person.offset{9} =  [120; 0;  0];
person.offset{10} = [120; 0;  0];

person.offset{11} = [0; 350;  0]; % pelvis
person.offset{12} = [-200; 0;  0];
person.offset{13} = [-200; 0;  0];
person.offset{14} = [-150; 0;  0];
person.offset{15} = [200; 0;  0];
person.offset{16} = [200; 0;  0];
person.offset{17} = [150; 0;  0];

figure(4); clf; hold on
set(gcf,'color','white')

person_generate

axis equal
view(153,23)
axis off
zoom(4)

view([-90 -90])

pbaspect([2 1 4])
figuresize(8,12,'centimeters')
matlabfrag('fig/hatze-bottom','renderer','opengl')
