function [calcs]= foot(O14,i_m,lr,left_foot_ankle_length,left_foot_toes_length,left_foot_heel_length,left_foot_upper_height,left_foot_lower_height,left_foot_length)

%% Foot

N = 100; % number of disks
indf = 1:N;

%% Densities

gamma_i = 1480*(1-(0.0001)*(indf.^2)*(1-1100/1480));
 
%% Measurements

a=left_foot_ankle_length;
b=left_foot_toes_length;
c=left_foot_heel_length;
h2=left_foot_upper_height;
h1= left_foot_lower_height;
l=left_foot_length;

%% Calculations

l_i=a+(l-a)*indf./N;
b_i=a+(b-a)*indf./N;
c_i=a+(c-a)*indf./N;
v_i = (h2/N)*(b_i.*l_i/2+c_i.*l_i/2);

m_i=gamma_i.*v_i;
m_11=990*l*h1*(5*b+c)/36;
m_12=990*l*h1*(b+5*c)/36;
m_13=1100*l*h1*(b+c)/4;

m=m_11+m_12+m_13+sum(m_i);
v=m_11/990+m_12/990+m_13/1100+sum(v_i);
 
% Mass centroid:
xc = 0;
yc = ((m_11+m_12)*(-h2-3*h1/4)...
    +m_13*(-h2-h1/4) ...
    +sum(m_i.*(-indf.*h2/N)))/m;
zc = (m_11*l*(-2/3+(7*b+2*c)/(9*(5*b+c)))...
    +m_12*l*(b+8*c)/(9*(b+5*c))+m_13*l*(-2/3+(b+2*c)/(3*(b+c)))...
    +sum(m_i.*(-a/2-(2*l/3-(2*l/3-a/2)*indf./100+l_i.*(b_i+2*c_i)/(3*(b_i+c_i))))))/m;

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

 
disp('-------------------------')
if lr == 'l'
  disp('Left foot section')
elseif lr == 'r'
  disp('Right foot section')
end
disp('-------------------------')
fprintf('Mass:     %2.3f kg\n',sum(m))
fprintf('Volume:   %1.4f m^3\n',sum(v))
fprintf('Centroid: [ %2.0f , %2.0f , %2.0f ] mm\n',1000*xc,1000*yc,1000*zc)
%fprintf('Moments of inertia: [ %2.3f , %2.3f , %2.3f ] kg.m^2\n',Ip_x,Ip_y,Ip_z)

calcs = [sum(m),sum(v),xc,yc,zc];% Ip_x,Ip_y,Ip_z];

%% Plot
 
opt  = {'opacity',0.1,'edgeopacity',0.1};
optl = {'opacity',0.2,'edgeopacity',0.1};
 
for ii = indf
  ph = -ii*h2/100; % plate height
  plot_trapzoidal_plate(O14+[0;ii/100*h2;ph],c_i(ii),b_i(ii),l_i(ii),h2/100)
end

end
 
