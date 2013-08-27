function person = segment_foot(person,S)

P = person.origin{S} + person.offset{S};
i_m = person.sex;

%% Foot

N = 100; % number of disks
ind = 1:N;
gamma = person.density.foot(ind);
 
%% Measurements

a  = person.meas{S}.all(1);
b  = person.meas{S}.all(2);
c  = person.meas{S}.all(3);
h1 = person.meas{S}.all(4);
h2 = person.meas{S}.all(5) - person.meas{S}.all(4);
l  = person.meas{S}.all(6);

%% Calculations

l_i=a+(l-a)*ind./N;
b_i=a+(b-a)*ind./N;
c_i=a+(c-a)*ind./N;

v = (h2/N)*(b_i.*l_i/2+c_i.*l_i/2);
v_11=990*l*h1*(5*b+c)/36;
v_12=990*l*h1*(b+5*c)/36;
v_13=1100*l*h1*(b+c)/4;
volume = v_11+v_12+v_13+sum(v);

m=gamma.*v;
m_11=990*v_11;
m_12=990*v_12;
m_13=1100*v_13;
mass = m_11+m_12+m_13+sum(m);
 
% Mass centroid:
xc = 0;
yc = ((m_11+m_12)*(-h2-3*h1/4)...
    +m_13*(-h2-h1/4) ...
    +sum(m.*(-ind.*h2/N)))/mass;
zc = (m_11*l*(-2/3+(7*b+2*c)/(9*(5*b+c)))...
    +m_12*l*(b+8*c)/(9*(b+5*c))+m_13*l*(-2/3+(b+2*c)/(3*(b+c)))...
    +sum(m.*(-a/2-(2*l/3-(2*l/3-a/2)*ind/N+l_i.*(b_i+2*c_i)/(3*(b_i+c_i))))))/mass;

% Moments of inertia:

%Ix_i=m_i.*l^2*(b_i^2+4*b_i.*c_i+c_i^2)/(18*(b_i+c_i)^2);
%Iz_i=m_i.*(b_i.^2+c_1.^2)/24;
%Iy_i=Ix_i+Iz_i;

% % principal moments of inertia; 
% Ip_x = m_11*((l/3)^2*(b^2+4*b*(2*b+c)/3+((2*b+c)/3)^2)/(18*(b+(2*b+c)/3)^2)+(h2+3*h1/4+yc)^2+(2*l/3-l*(7*b+2*c)/(9*(5*b+c))+zc)^2)...
%     +m_12*((l/3)^2*((b/3+2*c/3)^2+4*c*(b/3+2*c/3)+c^2)/(18*(b/3+5*c/3)^2)+(h2+3*h1/4+yc)^2+(l*(b+8*c)/(9*(b+5*c))-zc)^2)...
%     +m_13*(l^2*(b^2+4*b*c+c^2)/(18*(b+c)^2)+(h2+h1/4+yc)^2+(2*l/3-l*(b+c)+zc)^2)...
%     +sum(Ix_i+m_i.*((indf.*h2/N+yc).^2+a/2+(2*l/3-a/2)*inf./N-l_i.*(b_i+2*c_i)/(3*(b_i+c_i))+zc)^2);
% Iy1= m_11*((l/3)^2*(b^2+4*b*(2*b+c)/3+((2*b+c)/3)^2)/(18*(b+(2*b+c)/3)^2)+(b^2+(2*b/3+c/3)^2)/24+(2*l/3-l*(7*b+2*c)/(9*(5*b+c))+zc)^2)...
%      +m_12*((l/3)^2*((b/3+2*c/3)^2+4*c*(b/3+2*c/3)+c^2)/(18*(b/3+5*c/3)^2)+(c^2+(b/3+2*c/3)^2)/24+(l*(b+8*c)/(9*(b+5*c))-zc)^2)...
%      +m_13*(l^2*(b^2+4*b*c+c^2)/(18*(b+c)^2)+(b^2+c^2)/24+(2*l/3-l*(b+c)+zc)^2)...
%      +sum(Iy_i+m_i.*(a/2+(2*l/3-a/2)*inf./N-l_i.*(b_i+2*c_i)/(3*(b_i+c_i))+zc)^2);
% Iz1= m_11*(((b^2+(2*b+c)/3)^2)/24+(h2+3*h1/4+yc)^2)...
%     +m_12*(((b/3+2*c/3)^2+c^2)/24+(h2+3*h1/4+yc)^2)...
%     +m_13*((b^2+c^2)/24+(h2+h1/4+yc)^2)...
%     +sum(Iz_i+m_i.*(inf.*h2/N+yc)^2);
% Iy1z1= m_11*(-h2-3*h1/4-yc)*(-2*l/3+l*(7*b+2*c)/(9*(5*b+c))-zc)...
%     +m_12*(-h2-3*h1/4-yz)*(l*(b+8*c)/(9*(b+5*c))-zc)...
%     +m_13*(-h2-h1/4-c)*(-2*l/3+l*(b+2*c)/(3*(b+c))-zc)...
%     +sum(m_i.*(-indf.*h2/N-yc)*(-a/2-(2*l/3-a/2)*indf./N+l_i.*(b_i+2*c_i)/(3*(b_i+c_i)-zc)));
% 
% Ip_y = (Iy1+Iz1)/2+((Iy1-Iz1)^2/4+Iy1z1^2)^0.5;
% Ip_z = (Iy1+Iz1)/2-((Iy1-Iz1)^2/4+Iy1z1^2)^0.5;


person.segment(S).mass = mass;
person.segment(S).volume = volume;
person.segment(S).centroid = [xc; yc; zc];
%person.segment(S).Minertia = [Ip_x,Ip_y,Ip_z];
person.segment(S).theta = 30; % needs to be calculated

%% Plot

if person.plot
  
  opt1  = {'opacity',person.opacity{S}(1),'edgeopacity',0,'colour',person.color{S}};
  opt  = {'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2),'colour',person.color{S}};
  
  for ii = ind
    ph = -ii*h2/100; % plate height
    s = -a/2-ii/100*(2/3*l-a/2);
    plot_trapzoidal_plate(P+[0;-s-l_i(ii)/2;ph],b_i(ii),c_i(ii),l_i(ii),h2/100,opt1{:})
  end
  
  plot_trapzoidal_plate(P+[0;-s-l/2;-h2],b,c,l,-h1/2,opt{:})
  plot_trapzoidal_plate(P+[0;-s-l/2-l/3;-h2-h1],c+1/3*(b-c),c,l/3,h1/2,opt{:})
  plot_trapzoidal_plate(P+[0;-s-l/2+l/3;-h2-h1],b,b+1/3*(c-b),l/3,h1/2,opt{:})
  
end

end

