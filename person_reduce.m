%% Driver file to experiment with interpolating limbs
%
% NB that Nmeas interpolates down from Hatze's original measurements.
% The results of this interpolation is then used for interpolating to
% an arbitrary number (Ncalc) of larger values.

% Initialise

clear all
clc

pdata = @(person,s) [1000*person.segment(s).volume, person.segment(s).mass, person.segment(s).centroid(3), person.segment(s).Minertia(1), person.segment(s).Minertia(2), person.segment(s).Minertia(3)];

%% Original data

person = person_generate('data','hatze_meas.txt');

calc = [];
calc(5,end+1,:)  = pdata(person,5);
calc(9,end+1,:)   = pdata(person,9);
calc(4,end+1,:)        = pdata(person,4);
calc(5,end+1,:)       = pdata(person,5);

disp(squeeze(calc))

%% Figure showing interpolation

willfig('interp'); clf; hold on
figuresize(15,8,'centimeters');

forearm_right = 9;

person.segment(forearm_right).plot = true;

x = 250;

col = person.color{forearm_right};

person.offset{forearm_right} = [0; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 10; % SAME
person.segment(forearm_right).Ncalc = 10;
person.color{forearm_right} = [1 0 0];
person_generate(person)

person.offset{forearm_right} = [2*x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 10; % SAME
person.segment(forearm_right).Ncalc = 40;
person.color{forearm_right} = col;
person_generate(person)

person.offset{forearm_right} = [x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 5;
person.segment(forearm_right).Ncalc = 5;
person.color{forearm_right} = [1 0 0];
person_generate(person)

person.offset{forearm_right} = [3*x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 5;
person.segment(forearm_right).Ncalc = 40;
person.color{forearm_right} = col;
person_generate(person)

text(-0.19,0,0.1,'10 measurements','HorizontalAlignment','center')
text(-0.19+1*x/1000,0,0.1,'5 measurements','HorizontalAlignment','center')
text(-0.19+2*x/1000,0,0.1,'Interp. from 10','HorizontalAlignment','center')
text(-0.19+3*x/1000,0,0.1,'Interp. from 5','HorizontalAlignment','center')

view([0 0])
axis tight equal off

%matlabfrag('fig/interp-slice','renderer','opengl')


%% interpolation

clear all
clc

pdata = @(person,s) [1000*person.segment(s).volume, person.segment(s).mass, person.segment(s).centroid(3), person.segment(s).Minertia(1), person.segment(s).Minertia(2), person.segment(s).Minertia(3)];

person = person_generate('data','hatze_meas.txt');
forearm_right = 9;

person.segment(forearm_right).plot = false;

calcs = 3:10;
NC = numel(calcs);

person.segment(forearm_right).Nmeas = 10;
person.segment(forearm_right).Ncalc = 10;
person = person_generate(person);
forearm_R_standard = pdata(person,forearm_right);

for cc = 1:length(calcs)
  
  person.segment(forearm_right).Nmeas = calcs(cc);
  person.segment(forearm_right).Ncalc = calcs(cc);
  person = person_generate(person);
  
  forearm_R_normal(cc,:)  = abs(100*(pdata(person,forearm_right)-forearm_R_standard)./forearm_R_standard);

  person.segment(forearm_right).Nmeas = calcs(cc);
  person.segment(forearm_right).Ncalc = 10;
  person = person_generate(person);
  
  forearm_R_interp(cc,:)  = abs(100*(pdata(person,forearm_right)-forearm_R_standard)./forearm_R_standard);

end

willfig('interp results','small'); clf; hold on

plot(1:NC,forearm_R_normal(:,1),'-v')
plot(1:NC,forearm_R_normal(:,2),'-s')
plot(1:NC,forearm_R_normal(:,3),'-o')
plot(1:NC,forearm_R_normal(:,4),'--v')
plot(1:NC,forearm_R_normal(:,5),'--s')
plot(1:NC,forearm_R_normal(:,6),'--o')
colourplot
%plot(1:NC,sqrt(mean(forearm_R_normal.^2,2)),'k-','linewidth',2)
xlim([0 NC+1])
ylim([-0.5 4.5])
set(gca,'xtick',1:NC,'xticklabel',calcs,'ytick',[0:4],'yticklabel',[0:4])
xlabel('Number of slices')
ylabel('Percent error')

hleg = legend('$V$','$M$','$\bar z$','$\bar I_x^P$','$\bar I_y^P$','$\bar I_z^P$','location','northeast');
set(hleg,'box','off','interpreter','latex','fontsize',12)
%matlabfrag('fig/err-decimate')


willfig('interp2 results','small'); clf; hold on

plot(1:NC,forearm_R_interp(:,1),'-v')
plot(1:NC,forearm_R_interp(:,2),'-s')
plot(1:NC,forearm_R_interp(:,3),'-o')
plot(1:NC,forearm_R_interp(:,4),'--v')
plot(1:NC,forearm_R_interp(:,5),'--s')
plot(1:NC,forearm_R_interp(:,6),'--o')
colourplot
%plot(1:NC,sqrt(mean(forearm_R_interp.^2,2)),'k-','linewidth',2)
xlim([0 NC+1])
ylim([-0.5 4.5])
set(gca,'xtick',1:NC,'xticklabel',calcs,'ytick',[0:4],'yticklabel',[0:4])
xlabel('Number of slices before interpolation')
ylabel('Percent error')

hleg = legend('$V$','$M$','$\bar z$','$\bar I_x^P$','$\bar I_y^P$','$\bar I_z^P$','location','northeast');
set(hleg,'box','off','interpreter','latex','fontsize',12)
%matlabfrag('fig/err-interp')


