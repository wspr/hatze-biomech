%% Plot shoulder

clear all
figure(1); clf; hold on

person = person_generate('data','hatze_meas.txt');
person.segment(3).plot = true;
person.segment(3).opacity = [0.4 1];
person.segment(3).colour = [0 0.7 0.3];
person_generate(person)

axis equal
view([0 0])
pbaspect([1 1 1])
camlight
axis off

%%

figuresize(8,8,'centimeters')

view([0 0])
matlabfrag('fig/shoulder-front','renderer','opengl')

view([90 0])
matlabfrag('fig/shoulder-side','renderer','opengl')

view([0 90])
matlabfrag('fig/shoulder-top','renderer','opengl')