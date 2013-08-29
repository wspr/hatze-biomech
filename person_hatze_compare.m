%% HATZE COMPARE

clear all
clc

person_setup
person = person_data(person,'hatze_meas.txt');
person = person_generate(person);

%% Hatze's results to verify

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

%% Print results

for s = 1:person.N
  
  if ~isempty(person.segment(s).volume)
    disp('-------------------------')
    disp(person.segment(s).name)
    disp('-------------------------')
    fprintf('Volume:   %1.4f L\n',1000*person.segment(s).volume)
    fprintf('         (%1.4f)\n',person.segment(s).volume_hatze)
    fprintf('Mass:     %2.3f kg\n',person.segment(s).mass)
    fprintf('         (%2.3f)\n',person.segment(s).mass_hatze)
    if ~isempty(person.segment(s).centroid)
      fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*person.segment(s).centroid(1),1000*person.segment(s).centroid(2),1000*person.segment(s).centroid(3))
    end
    if ~isempty(person.segment(s).centroid_hatze)
      fprintf('         ([ %2.0f , %2.0f , %2.0f ])\n',1000*person.segment(s).centroid_hatze(1),1000*person.segment(s).centroid_hatze(2),1000*person.segment(s).centroid_hatze(3))
    end
    if ~isempty(person.segment(s).Minertia)
      fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] g.m^2\n',1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(3))
    end
    if ~isempty(person.segment(s).Minertia_hatze)
      fprintf('                   ([ %2.3f , %2.3f , %2.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(3))
    end
    
  end

end

