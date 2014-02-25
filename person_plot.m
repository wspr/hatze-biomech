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

person_generate(person,'plot',true);

person.plot_points([person.segment(3:6).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(7:10).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(11:14).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([11,15:17]).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([1,3,2]).origin], 'k.-', 'markersize', 20,'linewidth',2)

for ii = 1:person.N
  if ~isempty(person.segment(ii).centroid)
%    person.plot_points(person.segment(ii).origin+person.segment(ii).centroid, 'r.', 'markersize', 30)
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

person.plot_points([person.segment(3:6).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(7:10).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment(11:14).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([11,15:17]).origin], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.segment([1,3,2]).origin], 'k.-', 'markersize', 20,'linewidth',2)

for ii = 1:person.N
  if ~isempty(person.segment(ii).centroid)
    person.plot_points(person.segment(ii).origin+person.segment(ii).centroid, 'r.', 'markersize', 30)
  end
end

ind = 1:person.N;
ind([3 7]) = []; % repeated
for ii = ind
  plot_coord(person.segment(ii).origin,'index',[num2str(ii),''''],'rotate',person.segment(ii).Rglobal,'length',0.07);
end

%% Plot FRONT and SIDE views

% offsets for visualisations below

person.segment(1).offset =  [0; 0;  0];
person.segment(2).offset =  [0; 0; 50]/1000;
person.segment(3).offset =  [-100; 250; 100]/1000; % left shoulder
person.segment(4).offset =  [-50; 0; -140]/1000;
person.segment(5).offset =  [0; 0; -50]/1000;
person.segment(6).offset =  [0; 60; -50]/1000;
person.segment(7).offset =  [ 100; -250; 100]/1000; % right shoulder
person.segment(8).offset =  [50; 0; -140]/1000;
person.segment(9).offset =  [0; 0; -50]/1000;
person.segment(10).offset = [0; -60; -50]/1000;

person.segment(11).offset = [0; 0; -50]/1000; % pelvis
person.segment(12).offset = [0; 120; -140]/1000;
person.segment(13).offset = [0; 0; -50]/1000;
person.segment(14).offset = [0; 0; -50]/1000;
person.segment(15).offset = [0; -120; -140]/1000;
person.segment(16).offset = [0; 0; -50]/1000;
person.segment(17).offset = [0; 0; -50]/1000;

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



%% Plotting a single segment only

close all
figure(1); clf; hold on

person = person_generate('data','hatze_meas.txt');
person.segment(3).plot = true;
person.segment(3).opacity = [1 0];
person.segment(3).colour = [0 0.7 0.3];
person_generate(person)

axis equal
view([20 70])
pbaspect([1 1 1])
camlight