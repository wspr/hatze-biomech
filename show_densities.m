% densities

person = person_generate('data','hatze_meas_c_larkin.txt');

%%

willfig('thighhead'); clf; hold on
figuresize(13,6,'centimeters')

nu = linspace(0,0.2);
plot(nu,person.density.thigh_head(1,nu))
plot(nu,person.density.thigh_head(0,nu))
plot(nu,person.density.buttocks(1,nu))
plot(nu,person.density.buttocks(0,nu))

h = legend('Thigh head, M','Thigh head, F','Buttocks, M','Buttocks, F','location','eastoutside');
set(h,'box','off')
colourplot
xlabel('Subcutaneous fat indicator')
ylabel('Density')

matlabfrag('fig/thigh-head')

%%

willfig('thigh','small'); clf; hold on

for mf = [1 0]
nu = 0.05;

for ii = 1:10
  
  th(ii) = person.density.thigh{ii}(mf,nu);
  
end

d = ((1:10)-1)/9;

plot(d,th,'.-')
plot(-d(2),person.density.thigh_head(mf,nu),'o')

end

colourplot(2,[1 3 2 4])

xlim([-d(3) d(end)+d(2)])
xlabel('Thigh length')
ylabel('Density')
title('$\nu=0.05$')

matlabfrag('fig/thigh','epspad',[0 10 0 10])


%%

willfig('foot','small'); clf; hold on

n = 1:100;
ft = person.density.foot(n);
  
plot((n-1)/99,ft)

colourplot

xlabel('Foot height')
ylabel('Density')
title('Foot')

matlabfrag('fig/foot','epspad',[0 10 0 10])
