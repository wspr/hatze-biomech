%% Make a person

clear all
clc

person = person_generate('data','hatze_meas.txt');

%% Count measurements

count = 0;
fprintf('\n========================\nMeasurements per segment\n------------------------\n')
for ii = 1:person.N
    fprintf('%s: %i\n',person.segment(ii).name,numel(person.meas{ii}.all))
    %fprintf('%s & %i\\\\\n',person.segment(ii).name,numel(person.meas{ii}.all))
    count = count + numel(person.meas{ii}.all);
end
fprintf('-----------------------\nTotal measurements: %i\n=======================\n\n',count)

%% Plot isometric view of the model

person = person_generate('data','hatze_meas.txt');

figure(1); clf; hold on
set(gcf,'color','white')

axis equal
view(153,23)
axis off
zoom(2)

person = person_generate(person,'plot',true);
  
person.plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',2)

for ii = 1:person.N
  if ~isempty(person.segment(ii).centroid)
    person.plot_points(person.origin{ii}+person.segment(ii).centroid, 'r.', 'markersize', 30)
  end
end

pbaspect([2 1 4])

%figuresize(8,14,'centimeters')
%matlabfrag('fig/hatze-iso','renderer','opengl')

%% Plot skeleton view

figure(11); clf; hold on
set(gcf,'color','white')

axis equal
view(153,23)
axis off
pbaspect([2 1 4])
zoom(2)
  
person.plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',2)

for ii = 1:person.N
  if ~isempty(person.segment(ii).centroid)
    person.plot_points(person.origin{ii}+person.segment(ii).centroid, 'r.', 'markersize', 30)
  end
end

ind = 1:person.N;
ind([3 7]) = []; % repeated
for ii = ind
  plot_coord(person.origin{ii},'index',[num2str(ii),'''']);
end

%% Plot FRONT and SIDE views

% offsets for visualisations below

person.offset{1} =  [0; 0;  0];
person.offset{2} =  [0; 0; 50]/1000;
person.offset{3} =  [-100; 250; 100]/1000; % left shoulder
person.offset{4} =  [-50; 0; -140]/1000;
person.offset{5} =  [0; 0; -50]/1000;
person.offset{6} =  [0; 60; -50]/1000;
person.offset{7} =  [ 100; -250; 100]/1000; % right shoulder
person.offset{8} =  [50; 0; -140]/1000;
person.offset{9} =  [0; 0; -50]/1000;
person.offset{10} = [0; -60; -50]/1000;

person.offset{11} = [0; 0; -50]/1000; % pelvis
person.offset{12} = [0; 120; -140]/1000;
person.offset{13} = [0; 0; -50]/1000;
person.offset{14} = [0; 0; -50]/1000;
person.offset{15} = [0; -120; -140]/1000;
person.offset{16} = [0; 0; -50]/1000;
person.offset{17} = [0; 0; -50]/1000;

figure(2); clf; hold on
set(gcf,'color','white')

person_generate(person,'plot',true);

axis equal
view(153,23)
axis off
zoom(4)

view([180 0])

pbaspect([1 4 4])

%figuresize(8,16,'centimeters')
%matlabfrag('fig/hatze-front','renderer','opengl')


figure(3); clf; hold on
set(gcf,'color','white')

person_generate(person,'plot',true);

axis equal
view(153,23)
axis off
zoom(4)

view([90 0])
pbaspect([1 4 4])

%figuresize(8,16,'centimeters')
%matlabfrag('fig/hatze-side','renderer','opengl')


