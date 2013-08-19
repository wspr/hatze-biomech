function person = hand(person,S)

P = person.origin{S} + person.offset{S};

r = person.meas{S}(1);
h = person.meas{S}(2);

a10 = person.meas{S-1}.a(10);
b10 = person.meas{S-1}.b(10);

hh = 0.378*h;

opt = {'colour',person.color{S},'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2)};

plot_rect_prism(P,[2*a10 2*b10],[h 2.52*b10],-hh,opt{:});
plot_cylinder_hollow(P+[h/2;0;-hh-r],r-h/4,r,h,[0 180],'N',10,'rotate',[0 90 180],opt{:})
if S == 6 % left
  s = 3*h/4;
else
  s = 0;
end
plot_cylinder_hollow(P+[h/2-s;0;-hh-r],r-h/4,r,h/4,[180 180+180/pi*h/r],'N',10,'rotate',[0 90 180],opt{:})