%% Make a person

clear all
clc

person_data

count = 0;
for ii = 1:person.N
    fprintf('Segment: %s\nMeasurements: %i\n',person.segment(ii).name,numel(person.meas{ii}.all))
%    fprintf('%s & %i\\\\\n',person.segment(ii).name,numel(person.meas{ii}.all))
    count = count + numel(person.meas{ii}.all);
end
%count

person.plot = false;
person_generate

disp('=============')
disp('=============')
disp('=============')
disp('=============')
disp('=============')
disp('=============')

for s = 1:person.N
  
if ~isempty(person.segment(s).mass)
disp('-------------------------')
disp(person.segment(s).name)
disp('-------------------------')
fprintf('Volume:   %1.4f L\n',1000*person.segment(s).volume)
fprintf('         (%1.4f)\n',person.segment(s).volume_hatze)
fprintf('Mass:     %2.3f kg\n',person.segment(s).mass)
fprintf('         (%2.3f)\n',person.segment(s).mass_hatze)
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*person.segment(s).centroid(1),1000*person.segment(s).centroid(2),1000*person.segment(s).centroid(3))
fprintf('         ([ %2.0f , %2.0f , %2.0f ])\n',1000*person.segment(s).centroid_hatze(1),1000*person.segment(s).centroid_hatze(2),1000*person.segment(s).centroid_hatze(3))
fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] g.m^2\n',1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(3))
fprintf('                   ([ %2.3f , %2.3f , %2.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(3))
end

end

%% Calculate & plot

person.plot = true;

figure(1); clf; hold on
set(gcf,'color','white')

axis equal
view(153,23)
axis off
zoom(2)

person_generate
  
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
