function person = segment_hand(person,S)

P = person.origin{S} + person.offset{S};
person.segment(S).Rglobal = person.segment(S-1).Rglobal*person.segment(S).Rlocal;

r = person.meas{S}.all(1);
h = person.meas{S}.all(2);

a10 = person.meas{S-1}.a(end);
b10 = person.meas{S-1}.b(end);

hh = 0.378*h;

opt = {'colour',person.color{S},'opacity',person.opacity{S}(1),'edgeopacity',person.opacity{S}(2)};

%%

% % Mass 
% m_p = gamma*b(10)*(0.378*h)*(2.172*a(10)+1.172*h);
% m_c = gamma*pi*h^2*(2*r-h/4)/8;
% m_T = gamma*h^3/16;
% m = m_p+m_c+m_T;
% 
% v=1110*m/gamma;
% 
% % Mass centroid
% xc1=0.4244*(r^3-(r-h/4)^3)/(r^2-(r-h/4)^2);
% zp1=r*b(10)*(0.378*h)^2*(0.7533*a(10)+0.7967*h)/m_p;
% 
% xc = (¡À)*(m_c*xc1-m_T(r-h/4))/m;
% yc= m_T*3*h/8/m;
% zc=-(m_p*zp1+m_c*(r+0.378*h)+m_T*(r/3+0.378*h))/m;
% 
% % Moments of inertia:
% z1=-(zc+0.378*h);
% x1=-(xc+b(10)*(1-0.688*(z+zc)/h));
% x2=b(10)*(1-0.688(z+zc)/h)-xc;
% 
% Ig_xc=m_c*(r^2+(r-h/4)^2+h^2/3)/4;
% Ig_xy=m_c*(r^2+(r-h/4)^2)/2-xc1^2;
% Ig_xz=I_xcg-m_c*xc1^2;
% 
% I_x1=1.333*r*b(10)*(a(10)^3*0.378*h)+(3*a(10)^2*c(2)+...
%     c(1)*a(10)^3)*(0.378*h)^2/2+3*(a(10)*c(2)^2+...
%     c(1)*a(10)^2*c(2))*(0.378*h)^3/3+(c(2)^3+...
%     3*c(1)*a(10)*c(2)^2)*(0.378*h)^4/4+...
%     c(1)*c(2)^2*(0.378*h)^5/5);
% I_y2=1.333*r*b(10)^3*(a(10)*(0.378*h)+...
%     (3*a(10)*c(1)+c(2)*(0.378*h)^2/2+...
%     3*(a(10)*c(1)^2*+c(1)*c(2))*(0.378*h)^3/3+...
%     (a(10)*c(1)^3+3*c(1)*c(2)^2)*(0.378*h)^4/4+...
%     c(1)^3*c(2)*(0.378*h)^5/5);
% Ip=r*b(10)*(0.378*h)^3*(0.385*a(10)+0.604*h)-m_p*zp1^2;
% 
% I_x = I_x1+Ip+m_p*(yc^2+(zp1+zc^2)+Ig_xc+m_c*(yc^2+(r+0.378*h+zc)^2+m_T*((3*h/8-yc)^2+(r/3+0.378*h+zc)^2);
% I_z=I_x1+I_y2+m_p*(xc^2+yc^2)+Ig_xc+m_c*(yc^2+(xc1-abs(xc))^2-xc1^2)+m_T*((abs(xc)+r-h/4)^4+(3*h/8-yc)^2);
% I_xz=(¡À)r*h*((r^4-(r-h/4)^4)((sin(h/(h/8-r)))^2/32-((r^3-(r-h/4)^3)3)*((xc/4)*sin(h/(h/8-r))+...
%     (zc+r+0.378*h)*(2.25-0.25*cos(h/(h/8-r)))+...
%     0.5*(r^2-(r-h/4)^2*xc*(zc+r+0.378*h)*(pi+0.25*h/(h/8-r)))¡Àm_p*(-xc)*(-zc-zp1);
% 
% % principal moments of inertia; 
% Ip_y =I_y2+m_p*(xc^2+(zp1+zc))^2+Ig_yc+m_c*(xc1-abs(xc))^2+(r+0.378*h+zc)^2+m_T*((abs(xc)+r-h/4)^2+(r/3+0.378*h+zc)^2);
% Ip_x =(I_x+I_z)/2+((I_x-I_z)^2/4+I_xz^2)^(1/2);
% Ip_z =(I_x+I_z)/2-((I_x-I_z)^2/4+I_xz^2)^(1/2);

person.segment(S).theta = 4; % needs to be calculated

%%
if person.plot || person.segment(S).plot
  
  plot_rect_prism(P,[2*a10 2*b10],[h 2.52*b10],-hh,opt{:});
  plot_cylinder_hollow(P+[h/2;0;-hh-r],r-h/4,r,h,[0 180],'N',10,'rotate',[0 90 180],opt{:})
  if S == 6 % left
    s = 3*h/4;
  else
    s = 0;
  end
  plot_cylinder_hollow(P+[h/2-s;0;-hh-r],r-h/4,r,h/4,[180 180+180/pi*h/r],'N',10,'rotate',[0 90 180],opt{:})
  
end

end