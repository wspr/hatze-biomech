
clear all
clc

pdata = @(person,s) [1000*person.segment(s).volume, person.segment(s).mass, person.segment(s).centroid(3), person.segment(s).Minertia(1), person.segment(s).Minertia(2), person.segment(s).Minertia(3)];

person_setup
person_data

%% Original data

person_generate

calc = [];
calc(forearm_left,end+1,:)  = pdata(person,forearm_left);
calc(forearm_right,end,:)   = pdata(person,forearm_right);
calc(arm_left,end,:)        = pdata(person,arm_left);
calc(arm_right,end,:)       = pdata(person,arm_right);

%% Figure showing interpolation

willfig('interp'); clf; hold on

person.segment(forearm_right).plot = true;

x = 150;

col = person.color{forearm_right};

person.offset{forearm_right} = [0; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 10; % SAME
person.segment(forearm_right).Ncalc = 10;
person.color{forearm_right} = [1 0 0];
person_generate

person.offset{forearm_right} = [x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 10; % SAME
person.segment(forearm_right).Ncalc = 40;
person.color{forearm_right} = col;
person_generate

person.offset{forearm_right} = [2*x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 5;
person.segment(forearm_right).Ncalc = 5;
person.color{forearm_right} = [1 0 0];
person_generate

person.offset{forearm_right} = [3*x; 0; 0]./1000;
person.segment(forearm_right).Nmeas = 5;
person.segment(forearm_right).Ncalc = 40;
person.color{forearm_right} = col;
person_generate

text(-0.19,0,0.1,'Original','HorizontalAlignment','center')
text(-0.19+1*x/1000,0,0.1,'Interp.','HorizontalAlignment','center')
text(-0.19+2*x/1000,0,0.1,'Decimated','HorizontalAlignment','center')
text(-0.19+3*x/1000,0,0.1,'Dec. + Interp.','HorizontalAlignment','center')

view([0 0])
axis tight equal off

matlabfrag('fig/interp-slice','renderer','opengl')


%% interpolation

clear all
clc

pdata = @(person,s) [1000*person.segment(s).volume, person.segment(s).mass, person.segment(s).centroid(3), person.segment(s).Minertia(1), person.segment(s).Minertia(2), person.segment(s).Minertia(3)];

person_setup
person_data

person.segment(forearm_right).plot = false;

calcs = 3:10;
NC = numel(calcs);

person.segment(forearm_right).Nmeas = 10;
person.segment(forearm_right).Ncalc = 10;
person_generate
forearm_R_standard = pdata(person,forearm_right);

for cc = 1:length(calcs)
  
  person.segment(forearm_right).Nmeas = calcs(cc);
  person.segment(forearm_right).Ncalc = calcs(cc);
  person_generate
  
  forearm_R_normal(cc,:)  = abs(100*(pdata(person,forearm_right)-forearm_R_standard)./forearm_R_standard);

  person.segment(forearm_right).Nmeas = calcs(cc);
  person.segment(forearm_right).Ncalc = 10;
  person_generate
  
  forearm_R_interp(cc,:)  = abs(100*(pdata(person,forearm_right)-forearm_R_standard)./forearm_R_standard);

end

willfig('interp results','small'); clf; hold on

plot(1:NC,forearm_R_normal(:,1),'-v')
plot(1:NC,forearm_R_normal(:,2),'-s')
plot(1:NC,forearm_R_normal(:,3),'-o')
plot(1:NC,forearm_R_normal(:,4),'--v')
plot(1:NC,forearm_R_normal(:,5),'--s')
plot(1:NC,forearm_R_normal(:,6),'--o')
xlim([0 NC+1])
ylim([-0.5 4.5])
set(gca,'xtick',1:NC,'xticklabel',calcs,'ytick',[0:4],'yticklabel',[0:4])
xlabel('Number of slices')
ylabel('Percent error')
colourplot

hleg = legend('$V$','$M$','$\bar z$','$\bar I_x^P$','$\bar I_y^P$','$\bar I_z^P$','location','northeast');
set(hleg,'box','off','interpreter','latex','fontsize',12)
matlabfrag('fig/err-decimate')


willfig('interp2 results','small'); clf; hold on

plot(1:NC,forearm_R_interp(:,1),'-v')
plot(1:NC,forearm_R_interp(:,2),'-s')
plot(1:NC,forearm_R_interp(:,3),'-o')
plot(1:NC,forearm_R_interp(:,4),'--v')
plot(1:NC,forearm_R_interp(:,5),'--s')
plot(1:NC,forearm_R_interp(:,6),'--o')
xlim([0 NC+1])
ylim([-0.5 4.5])
set(gca,'xtick',1:NC,'xticklabel',calcs,'ytick',[0:4],'yticklabel',[0:4])
xlabel('Number of slices before interpolation')
ylabel('Percent error')
colourplot

hleg = legend('$V$','$M$','$\bar z$','$\bar I_x^P$','$\bar I_y^P$','$\bar I_z^P$','location','northeast');
set(hleg,'box','off','interpreter','latex','fontsize',12)
matlabfrag('fig/err-interp')


