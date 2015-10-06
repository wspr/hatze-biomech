%% HATZE COMPARE

clear all
clc

person = person_generate('data','hatze_meas_lechner.txt');
roundn = @(x,n) round(10^(n-1)*x)*10^(-n+1);

%% Hatze's results to verify

% thorax
person.segment(1).volume_hatze = 0.019890;
person.segment(1).mass_hatze = 18.487;
person.segment(1).centroid_hatze = [0; 0.002; 0.203];
person.segment(1).Minertia_hatze = [0.331264 0.377450 0.150689];
person.segment(1).theta_hatze = 0.023;

% head-neck
person.segment(2).volume_hatze = 0.004669;
person.segment(2).mass_hatze   = 5.187;
person.segment(2).centroid_hatze = [0; 0; 0.137];
person.segment(2).Minertia_hatze = [0.032585 0.027423 0.018451];

% left shoulder
person.segment(3).volume_hatze = 0.001551;
person.segment(3).mass_hatze   = 1.628;
person.segment(3).centroid_hatze = [0; 0; 0.153];
person.segment(3).Minertia_hatze = [0.005280 0.003784 NaN];
person.segment(3).resting_angle_hatze = 0.196;

% right shoulder
person.segment(7).volume_hatze = 0.001977;
person.segment(7).mass_hatze   = 2.076;
person.segment(7).centroid_hatze = [0; 0; 0.158];
person.segment(7).Minertia_hatze = [0.007601 0.005615 0.008042];  %Ipz=NaN or 0 from p40 
person.segment(7).resting_angle_hatze = -0.141;

% left arm
person.segment(4).volume_hatze = 0.002142;
person.segment(4).mass_hatze = 2.320;
person.segment(4).centroid_hatze = [0; 0; -0.131];
person.segment(4).Minertia_hatze = [0.019671 0.019450 0.002687];

% right arm
person.segment(8).volume_hatze = 0.002180;
person.segment(8).mass_hatze = 2.362;
person.segment(8).centroid_hatze = [0; 0; -0.129];
person.segment(8).Minertia_hatze = [0.019906 0.019744 0.002791];

% left forearm
person.segment(5).volume_hatze = 0.001072;
person.segment(5).mass_hatze = 1.177;
person.segment(5).centroid_hatze = [0; 0; -0.112];
person.segment(5).Minertia_hatze = [0.007228 0.007425 0.000837];

% right forearm
person.segment(9).volume_hatze = 0.001229;
person.segment(9).mass_hatze = 1.343;
person.segment(9).centroid_hatze = [0; 0; -0.114];
person.segment(9).Minertia_hatze = [0.008281 0.008608 0.001083];

% left hand
person.segment(6).volume_hatze = 0.000488;
person.segment(6).mass_hatze = 0.542;
person.segment(6).centroid_hatze = [-0.012; 0.003; -0.061];
person.segment(6).Minertia_hatze = [0.000578 0.001149 0.001316];
person.segment(6).theta_hatze = -1.198;

% right hand
person.segment(10).volume_hatze = 0.000477;
person.segment(10).mass_hatze = 0.529;
person.segment(10).centroid_hatze = [0.012; 0.004; -0.063];
person.segment(10).Minertia_hatze = [0.000562 0.001021 0.001236];
person.segment(10).theta_hatze = 1.198;

% pelvis
person.segment(11).volume_hatze = 0.009044;
person.segment(11).mass_hatze = 9.479;
person.segment(11).centroid_hatze = [0; -0.005; -0.079];
person.segment(11).Minertia_hatze = [0.046296 0.101353 0.058274];
person.segment(11).theta_hatze = -1.352;

% left thigh
person.segment(12).volume_hatze = 0.008360;
person.segment(12).mass_hatze = 8.938;
person.segment(12).centroid_hatze = [0; 0; -0.217];
person.segment(12).Minertia_hatze = [0.150227 0.146865 0.028530];

% right thigh
person.segment(15).volume_hatze = 0.008344;
person.segment(15).mass_hatze = 8.915;
person.segment(15).centroid_hatze = [0; 0; -0.208];
person.segment(15).Minertia_hatze = [0.141516 0.140247 0.029110];

% left leg
person.segment(13).volume_hatze = 0.003671;
person.segment(13).mass_hatze = 3.997;
person.segment(13).centroid_hatze = [0; 0; -0.186];
person.segment(13).Minertia_hatze = [0.061723 0.061303 0.005846];

% right leg
person.segment(16).volume_hatze = 0.003752;
person.segment(16).mass_hatze = 4.089;
person.segment(16).centroid_hatze = [0; 0; -0.194];
person.segment(16).Minertia_hatze = [0.068201 0.067841 0.005836];

% left foot
person.segment(14).volume_hatze = 0.000968;
person.segment(14).mass_hatze = 1.098;
person.segment(14).centroid_hatze = [0; -0.040; -0.039];
person.segment(14).Minertia_hatze = [0.004711 0.004958 0.000947];
person.segment(14).theta_hatze = -0.063;
person.segment(14).contact_heel_hatze = [0; -0.78; 0.052];
person.segment(14).contact_toe_hatze = [0; -0.067; -0.135];

% right foot
person.segment(17).volume_hatze = 0.000979;
person.segment(17).mass_hatze = 1.109;
person.segment(17).centroid_hatze = [0; -0.038; -0.038];
person.segment(17).Minertia_hatze = [0.004700 0.004951 0.000967];
person.segment(17).theta_hatze = -0.066;
person.segment(17).contact_heel_hatze = [0; -0.078; 0.051];
person.segment(17).contact_toe_hatze = [0; -0.066; -0.125];

%% Print results


for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    disp('-------------------------')
    disp([num2str(s),': ',person.segment(s).name])
    disp('-------------------------')
    fprintf('Volume:   %1.3f L\n',roundn(1000*person.segment(s).volume,3))
    fprintf('         (%1.3f)\n',1000*person.segment(s).volume_hatze)
    fprintf('Mass:     %2.3f kg\n',roundn(person.segment(s).mass,3))
    fprintf('         (%2.3f)\n',person.segment(s).mass_hatze)
    if ~isempty(person.segment(s).centroid)
      fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*person.segment(s).centroid(1),1000*person.segment(s).centroid(2),1000*person.segment(s).centroid(3))
    end
    if ~isempty(person.segment(s).centroid_hatze)
      fprintf('         ([ %2.0f , %2.0f , %2.0f ])\n',1000*person.segment(s).centroid_hatze(1),1000*person.segment(s).centroid_hatze(2),1000*person.segment(s).centroid_hatze(3))
    end
    if ~isempty(person.segment(s).Minertia)
      fprintf('Moments of inertia: [ %3.3f , %3.3f , %3.3f ] g.m^2\n',roundn(1000*person.segment(s).Minertia(1),3),roundn(1000*person.segment(s).Minertia(2),3),roundn(1000*person.segment(s).Minertia(3),3))
    end
    if ~isempty(person.segment(s).Minertia_hatze)
      fprintf('                   ([ %3.3f , %3.3f , %3.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(2),1000*person.segment(s).Minertia_hatze(3))
    end
    if ~isempty(person.segment(s).theta)
      fprintf('Theta: %2.3f°\n',person.segment(s).theta*180/pi)
    end
    if ~isempty(person.segment(s).theta_hatze)
      fprintf('      (%2.3f°)\n',person.segment(s).theta_hatze*180/pi)
    end

  end

end


%% Print results

clc

thresh = 0.01;

fprintf('\n\n===============================\n')
fprintf('=== ERRORS OF MORE THAN %i%% ===\n',round(100*thresh))
fprintf('===============================\n\n\n')


errf = @(a,b) (a-b)./b;
errfn = @(a,b,n) errf(roundn(a,n),roundn(b,n));
errsn = @(s,p,n) errfn(person.segment(s).(p),person.segment(s).([p,'_hatze']),n);

for s = 1:person.N

  print_err = true;
  
  if ~isempty(person.segment(s).volume)
    err = errsn(s,'volume',6);
    person.segment(s).volume_err = err;
    if abs(err) > thresh
      if print_err; print_err = false; disp(['--- ',person.segment(s).name,' ---']), end
      fprintf('Volume:   %1.4f L\n',1e3*roundn(person.segment(s).volume,6))
      fprintf('         (%1.4f)\n', 1e3*roundn(person.segment(s).volume_hatze,6))
      fprintf('Error:   %1.2f%% compared to Hatze\n' ,100*err);
    end
    
    err = errsn(s,'mass',3);
    person.segment(s).mass_err = err;
    if abs(err) > thresh      
      if print_err; print_err = false; disp(['--- ',person.segment(s).name,' ---']), end
      fprintf('Mass:     %2.3f kg\n',roundn(person.segment(s).mass,3))
      fprintf('         (%2.3f)\n',  roundn(person.segment(s).mass_hatze,3))
      fprintf('Error:   %1.2f%% compared to Hatze\n' ,100*err)
    end
    
    err = errsn(s,'centroid',3);
    err(isnan(err)) = 0;
    person.segment(s).centroid_err = err;
    if any(abs(err) > thresh)
      if print_err; print_err = false; disp(['--- ',person.segment(s).name,' ---']), end
      fprintf('Centroid: [ %2.1f , %2.1f , %2.1f ] mm\n',1000*roundn(person.segment(s).centroid(1),3),1000*roundn(person.segment(s).centroid(2),3),1000*roundn(person.segment(s).centroid(3),3))
      fprintf('         ([ %2.1f , %2.1f , %2.1f ])\n',1000*person.segment(s).centroid_hatze(1),1000*person.segment(s).centroid_hatze(2),1000*person.segment(s).centroid_hatze(3))
      fprintf('Error: [ %2.2f%% , %2.2f%% , %2.2f%% ] compared to Hatze\n',100*err(1),100*err(2),100*err(3))
    end
    
    err = errsn(s,'Minertia',6);
    err(isnan(err)) = 0;
    person.segment(s).Minertia_err = err;
    if any(abs(err) > thresh)
      if print_err; print_err = false; disp(['--- ',person.segment(s).name,' ---']), end
      fprintf('Moments of inertia: [ %3.3f , %3.3f , %3.3f ] g.m^2\n',1000*roundn(person.segment(s).Minertia(1),6),1000*roundn(person.segment(s).Minertia(2),6),1000*roundn(person.segment(s).Minertia(3),6))
      fprintf('                   ([ %3.3f , %3.3f , %3.3f ])\n',1000*person.segment(s).Minertia_hatze(1),1000*person.segment(s).Minertia_hatze(2),1000*person.segment(s).Minertia_hatze(3))
      fprintf('Error: [ %2.2f%% , %2.2f%% , %2.2f%% ] compared to Hatze\n',100*err(1),100*err(2),100*err(3))
    end
    
    if ~isempty(person.segment(s).theta)
      err = errsn(s,'theta',6);
      err(isnan(err)) = 0;
      person.segment(s).theta_err = err;
      if abs(err) > thresh
        if print_err; print_err = false; disp(['--- ',person.segment(s).name,' ---']), end
        fprintf('Theta: %2.3f°\n',person.segment(s).theta*180/pi)
        fprintf('      (%2.3f°)\n',person.segment(s).theta_hatze*180/pi)
        fprintf('Error:   %1.2f%% compared to Hatze\n' ,100*err)
      end
    end
    
    if ~print_err; fprintf('\n\n\n'); end
  end

end

%% LaTeX print all

disp('\begin{tabular}{llccc}')
disp('\toprule')
disp('Segment & Parameter & Calculated & Given & Error (\%) \\')
  
for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    disp('\midrule')
    fprintf('%s & ',person.segment(s).name)
    fprintf('Volume (L) & ')
    fprintf('%1.3f & ',roundn(1000*person.segment(s).volume,3))
    fprintf('%1.3f & ',1000*person.segment(s).volume_hatze)
    fprintf('%1.1f \\\\ \n',100*person.segment(s).volume_err)
    fprintf('& Mass (kg) & ')
    fprintf('%2.3f & ',roundn(person.segment(s).mass,3))
    fprintf('%2.3f & ',person.segment(s).mass_hatze)
    fprintf('%1.1f \\\\ \n',100*person.segment(s).mass_err)
    for ii = 1:3
      if person.segment(s).centroid(ii) ~= 0
        fprintf('& Centroid, $%s$ (mm) & ',char('w'+ii))
        fprintf('%3.0f & ',1000*roundn(person.segment(s).centroid(ii),3))
        fprintf('%3.0f & ',1000*roundn(person.segment(s).centroid_hatze(ii),3))
        fprintf('%1.1f \\\\ \n',100*person.segment(s).centroid_err(ii))
      end
    end
    for ii = 1:3
      if person.segment(s).Minertia(ii) ~= 0
        fprintf('& Moment of inertia, $%s$ (\\si{g.m^2}) & ',char('w'+ii))
        fprintf('%3.3f & ',roundn(1000*person.segment(s).Minertia(ii),3))
        fprintf('%3.3f & ',roundn(1000*person.segment(s).Minertia_hatze(ii),3))
        fprintf('%1.1f \\\\ \n',100*person.segment(s).Minertia_err(ii))
      end
    end
    if ~isempty(person.segment(s).theta)
      fprintf('& Principle angle & %2.3f\\textdegree & ',person.segment(s).theta*180/pi)
      fprintf('%2.3f\\textdegree & ',person.segment(s).theta_hatze*180/pi)
      fprintf('%1.1f \\\\ \n ',person.segment(s).theta_err*100)
    end

  end

end

disp('\bottomrule')
disp('\end{tabular}')



%% LaTeX print condensed
%
% Needs the booktabs and siunitx packages.

disp('\def\colsp{\hspace{0.6em}}')
disp('\begin{tabular}{@{}lr@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r @{}}')
disp('\toprule')
disp('& \multicolumn{3}{c}{$M$, kg} & \multicolumn{3}{c}{$\bar z$, mm} & \multicolumn{3}{c}{$I_{x}$, \si{g.m^2}} & \multicolumn{3}{c}{$I_{y}$, \si{g.m^2}} & \multicolumn{3}{c}{$I_{z}$, \si{g.m^2}} \\')
disp('\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}\cmidrule(lr){11-13}\cmidrule(lr){14-16}')
disp('Segment & Calc & Hatze & Err, \% &  Calc & Hatze & Err, \% &  Calc & Hatze & Err, \% &  Calc & Hatze & Err, \% &  Calc & Hatze & Err, \%  \\')
disp('\midrule')

for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    fprintf('%s ',person.segment(s).name)
    fprintf('& \\num{%2.3f} ',roundn(person.segment(s).mass,4))
    fprintf('& \\num{%2.3f} ',person.segment(s).mass_hatze)
    fprintf('& \\num{%1.2f} ',100*person.segment(s).mass_err)
    for ii = 3
        fprintf('& \\num{%3.0f} ',1000*roundn(person.segment(s).centroid(ii),4))
        fprintf('& \\num{%3.0f} ',1000*roundn(person.segment(s).centroid_hatze(ii),4))
        fprintf('& \\num{%1.2f} ',100*person.segment(s).centroid_err(ii))

    end
    for ii = 1:3
        fprintf('& \\num{%3.3f} ',roundn(1000*person.segment(s).Minertia(ii),4))
        if isnan(person.segment(s).Minertia_hatze(ii))
          fprintf('& --- ')
        else
          fprintf('& \\num{%3.3f} ',roundn(1000*person.segment(s).Minertia_hatze(ii),4))
        end
        fprintf('& \\num{%1.2f} ',100*person.segment(s).Minertia_err(ii))
    end
    fprintf('\\\\\n')

  end

end

disp('\bottomrule')
disp('\end{tabular}')


%% LaTeX print condensed
%
% Needs the booktabs and siunitx packages.

clc

disp('\def\colsp{\hspace{0.6em}}')
disp('\begin{tabular}{@{}lr@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r r@{\colsp}r@{\colsp}r @{}}')
disp('\toprule')
disp('& \multicolumn{3}{c}{$M$, kg} & \multicolumn{3}{c}{$\bar z$, mm} & \multicolumn{3}{c}{$I_{x}$, \si{g.m^2}} \\')
disp('\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}')
disp('Segment & Calc & Hatze & Err, \% &  Calc & Hatze & Err, \% &  Calc & Hatze & Err, \%  \\')
disp('\midrule')

for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    fprintf('%s ',person.segment(s).name)
    fprintf('& \\num{%2.3f} ',roundn(person.segment(s).mass,4))
    fprintf('& \\num{%2.3f} ',person.segment(s).mass_hatze)
    fprintf('& \\num{%1.2f} ',100*person.segment(s).mass_err)
    for ii = 3
        fprintf('& \\num{%3.0f} ',1000*roundn(person.segment(s).centroid(ii),4))
        fprintf('& \\num{%3.0f} ',1000*roundn(person.segment(s).centroid_hatze(ii),4))
        fprintf('& \\num{%1.2f} ',100*person.segment(s).centroid_err(ii))

    end
    for ii = 1
        fprintf('& \\num{%3.3f} ',roundn(1000*person.segment(s).Minertia(ii),4))
        if isnan(person.segment(s).Minertia_hatze(ii))
          fprintf('& --- ')
        else
          fprintf('& \\num{%3.3f} ',roundn(1000*person.segment(s).Minertia_hatze(ii),4))
        end
        fprintf('& \\num{%1.2f} ',100*person.segment(s).Minertia_err(ii))
    end
    fprintf('\\\\\n')

  end

end

disp('\bottomrule')
disp('\end{tabular}')

%% LaTeX print condensed just mass
%
% Needs the booktabs and siunitx packages.

disp('\def\colsp{\hspace{0.6em}}')
disp('\begin{tabular}{@{}lr@{\colsp}r@{\colsp}r @{}}')
disp('\toprule')
disp('& \multicolumn{3}{c}{Mass, kg} \\')
disp('\cmidrule(l){2-4}')
disp('Segment & Calc & Hatze & Err, \%  \\')
disp('\midrule')

for s = 1:person.N

  if ~isempty(person.segment(s).volume)
    fprintf('%s ',person.segment(s).name)
    fprintf('& \\num{%2.3f} ',roundn(person.segment(s).mass,4))
    fprintf('& \\num{%2.3f} ',person.segment(s).mass_hatze)
    fprintf('& \\num{%1.2f} ',100*person.segment(s).mass_err)
    fprintf('\\\\\n')
  end

end

disp('\midrule')
fprintf('Total & \\num{%2.3f} & \\num{%2.3f} \\\\\n',sum([person.segment.mass]),sum([person.segment.mass_hatze]))
fprintf('Measured & \\multicolumn{2}{c}{\\num{%2.1f}} \\\\\n',person.weight)
disp('\bottomrule')
disp('\end{tabular}')

