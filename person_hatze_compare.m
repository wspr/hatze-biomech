%% HATZE COMPARE

clear all
clc

person = person_generate('data','hatze_meas.txt');

%% Hatze's results to verify

%data for R. Marga 31 year old female (Hatze 79)

% thorax
person.segment(1).volume_hatze = 14.537;
person.segment(1).mass_hatze = 12.950;
person.segment(1).centroid_hatze = [0; 0.003; 0.194];  
person.segment(1).Minertia_hatze = [0.206372 0.231513 0.082898];  %with respect to centroid
person.segment(1).theta_hatze = -0.001;  %angle between principle z-axis and original z-axis of segment

% head-neck
person.segment(2).volume_hatze = 3.595;
person.segment(2).mass_hatze   = 3.993;
person.segment(2).centroid_hatze = [0; 0; 0.132];
person.segment(2).Minertia_hatze = [0.020953 0.017553 0.01216];

% left shoulder
person.segment(3).volume_hatze = 1.11;
person.segment(3).mass_hatze   = 1.144;
person.segment(3).centroid_hatze = [0; 0; 0.140];
person.segment(3).Minertia_hatze = [0.002725 0.002172 NaN];
person.segment(3).theta_hatze = 0.199;  %resting inclination angle of segment z-axis to horizontal

% right shoulder
person.segment(7).volume_hatze = 1.146;
person.segment(7).mass_hatze   = 1.180;
person.segment(7).centroid_hatze = [0; 0; 0.141];
person.segment(7).Minertia_hatze = [0.002922 0.002289 NaN];
person.segment(7).theta_hatze = -0.182;

% left arm
person.segment(4).volume_hatze = 1.616;
person.segment(4).mass_hatze = 1.713;
person.segment(4).centroid_hatze = [0; 0; -0.131];
person.segment(4).Minertia_hatze = [0.013677 0.013435 0.001528];

% right arm
person.segment(8).volume_hatze = 1.505;
person.segment(8).mass_hatze = 1.595;
person.segment(8).centroid_hatze = [0; 0; -0.129];
person.segment(8).Minertia_hatze = [0.011917 0.011937 0.001325];

% left forearm
person.segment(5).volume_hatze = 0.835;
person.segment(5).mass_hatze = 0.896;
person.segment(5).centroid_hatze = [0; 0; -0.108];
person.segment(5).Minertia_hatze = [0.004596 0.004712 0.000523];

% right forearm
person.segment(9).volume_hatze = 0.809;
person.segment(9).mass_hatze = 0.869;
person.segment(9).centroid_hatze = [0; 0; -0.112];
person.segment(9).Minertia_hatze = [0.004765 0.004855 0.000472];

% left hand
person.segment(6).volume_hatze = 0.285;
person.segment(6).mass_hatze = 0.317;
person.segment(6).centroid_hatze = [-0.01; 0.003; -0.049];
person.segment(6).Minertia_hatze = [0.000231 0.000448 0.000521];
person.segment(6).theta_hatze = -1.185;

% right hand
person.segment(10).volume_hatze = 0.288;
person.segment(10).mass_hatze = 0.319;
person.segment(10).centroid_hatze = [0.01; 0.003; -0.51];
person.segment(10).Minertia_hatze = [0.000240 0.000491 0.000551];
person.segment(10).theta_hatze = -1.195;

% pelvis
person.segment(11).volume_hatze = 11.208;
person.segment(11).mass_hatze = 11.304;
person.segment(11).centroid_hatze = [0; 0.092; -0.050];
person.segment(11).Minertia_hatze = [0.068933 0.123422 0.08784];
person.segment(11).theta_hatze = -1.1;

% left thigh
person.segment(12).volume_hatze = 9.166;
person.segment(12).mass_hatze = 9.593;
person.segment(12).centroid_hatze = [0; 0; -0.194];
person.segment(12).Minertia_hatze = [0.156305 0.152523 0.035853];

% right thigh
person.segment(15).volume_hatze = 8.955;
person.segment(15).mass_hatze = 9.375;
person.segment(15).centroid_hatze = [0; 0; -0.194];
person.segment(15).Minertia_hatze = [0.149153 0.146883 0.034605];

% left leg
person.segment(13).volume_hatze = 3.310;
person.segment(13).mass_hatze = 3.602;
person.segment(13).centroid_hatze = [0; 0; -0.174];
person.segment(13).Minertia_hatze = [0.045315 0.044973 0.005106];

% right leg
person.segment(16).volume_hatze = 3.487;
person.segment(16).mass_hatze = 3.794;
person.segment(16).centroid_hatze = [0; 0; -0.183];
person.segment(16).Minertia_hatze = [0.050537 0.050114 0.005438];

% left foot
person.segment(14).volume_hatze = 0.842;
person.segment(14).mass_hatze = 0.961;
person.segment(14).centroid_hatze = [0; -0.034; -0.040];
person.segment(14).Minertia_hatze = [0.003610 0.003781 0.000763];
person.segment(14).theta_hatze = -0.08;
person.segment(14).contact_heel_hatze = [0; -0.76; 0.047];
person.segment(14).contact_toe_hatze = [0; -0.06; -0.155];

% right foot
person.segment(17).volume_hatze = 0.887;
person.segment(17).mass_hatze = 1.02;
person.segment(17).centroid_hatze = [0; -0.037; -0.036];
person.segment(17).Minertia_hatze = [0.003706 0.003847 0.000859];
person.segment(17).theta_hatze = -0.086;
person.segment(17).contact_heel_hatze = [0; -0.81; 0.046];
person.segment(17).contact_toe_hatze = [0; -0.065; -0.145];

%% Print results

for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    disp('-------------------------')
    disp([num2str(s),': ',person.segment(s).name])
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
      fprintf('Moments of inertia: [ %3.3f , %3.3f , %3.3f ] g.m^2\n',1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(2),1000*person.segment(s).Minertia(3))
    end
    if ~isempty(person.segment(s).Minertia_hatze)
      fprintf('                   ([ %3.3f , %3.3f , %3.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(2),1000*person.segment(s).Minertia_hatze(3))
    end
    if ~isempty(person.segment(s).theta)
      fprintf('Theta: %2.3f rad\n',person.segment(s).theta) %*180/pi)
    end
    if ~isempty(person.segment(s).theta_hatze)
      fprintf('      (%2.3f rad)\n',person.segment(s).theta_hatze) %*180/pi)
    end

  end

end


%% Print results

thresh = 0.1;

fprintf('\n\n===============================\n')
fprintf('=== ERRORS OF MORE THAN %i%% ===\n',round(100*thresh))
fprintf('===============================\n\n\n')


for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    if abs(1000*person.segment(s).volume-person.segment(s).volume_hatze)/person.segment(s).volume_hatze > thresh      
      disp(['--- ',person.segment(s).name,' ---'])
      fprintf('Volume:   %1.4f L\n',1000*person.segment(s).volume)
      fprintf('         (%1.4f)\n',person.segment(s).volume_hatze)
    end
    if abs(person.segment(s).mass-person.segment(s).mass_hatze)/person.segment(s).mass_hatze > thresh      
      disp(['--- ',person.segment(s).name,' ---'])
      fprintf('Mass:     %2.3f kg\n',person.segment(s).mass)
      fprintf('         (%2.3f)\n',person.segment(s).mass_hatze)
    end
    if any(abs(person.segment(s).centroid-person.segment(s).centroid_hatze)./person.segment(s).centroid_hatze > thresh)
      disp(['--- ',person.segment(s).name,' ---'])
      fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*person.segment(s).centroid(1),1000*person.segment(s).centroid(2),1000*person.segment(s).centroid(3))
      fprintf('         ([ %2.0f , %2.0f , %2.0f ])\n',1000*person.segment(s).centroid_hatze(1),1000*person.segment(s).centroid_hatze(2),1000*person.segment(s).centroid_hatze(3))
    end
    if any(abs(person.segment(s).Minertia-person.segment(s).Minertia_hatze)./person.segment(s).Minertia_hatze > thresh)
      disp(['--- ',person.segment(s).name,' ---'])
      fprintf('Moments of inertia: [ %3.3f , %3.3f , %3.3f ] g.m^2\n',1000*person.segment(s).Minertia(1),1000*person.segment(s).Minertia(2),1000*person.segment(s).Minertia(3))
      fprintf('                   ([ %3.3f , %3.3f , %3.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(2),1000*person.segment(s).Minertia_hatze(3))
    end
    if ~isempty(person.segment(s).theta)
      if abs( (person.segment(s).theta-person.segment(s).theta_hatze )/person.segment(s).theta_hatze) > thresh
        disp(['--- ',person.segment(s).name,' ---'])
        fprintf('Theta: %2.3f rad\n',person.segment(s).theta) %*180/pi)
        fprintf('      (%2.3f rad)\n',person.segment(s).theta_hatze) %*180/pi)
      end
    end
  end

end


