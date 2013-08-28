%% Driver file to experiment with movement

clear all
close all
clc


person_setup
person_data

person = person_generate(person);

person.q(4)  =  person.segment( 1).theta;
person.q(7)  = -person.segment( 1).theta;
person.q(10) = -person.segment( 1).theta;
person.q(18) =  person.segment( 6).theta;
person.q(19) = -person.segment( 1).theta;
person.q(27) =  person.segment(10).theta;
person.q(28) =  person.segment(11).theta - person.segment(1).theta;
person.q(31) = -person.segment(11).theta;
person.q(35) =  person.segment(14).theta + 90;
person.q(37) = -person.segment(11).theta;
person.q(41) =  person.segment(17).theta + 90;

person.q(4) = 0; % lumbar extension
person.q(6) = 0; % lumbar twist

person = person_generate(person,'plot',true);

person.plot_points([person.origin{3:6}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{7:10}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{11:14}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[11,15:17]}], 'k.-', 'markersize', 20,'linewidth',2)
person.plot_points([person.origin{[1,3,2]}], 'k.-', 'markersize', 20,'linewidth',2)

axis equal
view([140 20])
pbaspect([1 1 4])