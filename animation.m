%% animation.m
%Date:2020/11/12
%Code writter: Yi-Lun Hsu <r07921067@ntu.edu.tw> 
%Professor: Cheng-Wei Chen <cwchenee@ntu.edu.tw>
%All rights reserved.
%----------------------------------------------------------------------------------
% x [m]: x-position matrix of body
% z [m]: z-position matrix of body
% theta [rad]: angle matrix 
% q1x [m]: x-position matrix of footpad1
% q1z [m]: z-position matrix of footpad1
% q2x [m]: x-position matrix of footpad2
% q2z [m]: z-position matrix of footpad2
%--------------------------------------------------
function []=animation(x,z,theta,q1x,q1z,q2x,q2z)
N = length(x);
W = 1;         %(m)- Width of the body
H = 1;         %(m)- Height of the body

% Animation-----------------------------
%--------------RGB color-----------------
rgb =[0, 0.4470, 0.7410;...     %new_blue
      0.8500, 0.3250, 0.0980;...%orange
      0.9290, 0.6940, 0.1250;...%soil yellow
      0.4940, 0.1840, 0.5560;...%purple
      0.4660, 0.6740, 0.1880;...%lightly green
      0.3010, 0.7450, 0.9330;...%lightly blue
      0.6350, 0.0780, 0.1840];  %Brown
%----------------------------------------

%----------------------------------------
figure(2);
axis ([-10 10 0 10]); 
axis manual;

for n = 1:20:N
    Rc=[cos(theta(n)),-sin(theta(n));sin(theta(n)),cos(theta(n))];
    qb=[x(n);z(n)];
    
    r1_c = [-W/2;H/2];
    r2_c = [W/2;H/2]; 
    r3_c = [-W/2;-H/2];
    r4_c = [W/2;-H/2];
    
    r1 = Rc*r1_c;
    r2 = Rc*r2_c;
    r3 = Rc*r3_c;
    r4 = Rc*r4_c;
    
    q1 = qb + r1;
    q2 = qb + r2;
    q3 = qb + r3;
    q4 = qb + r4;
    
    plot(x(n),z(n),'ro');
    hold on;
    plot(q1x(n),q1z(n),'bo');
    plot(q2x(n),q2z(n),'bo');
    %--------------Lander body-------------
    plot([q1(1) q2(1)],[q1(2) q2(2)],'color',rgb(7,:));
    plot([q2(1) q4(1)],[q2(2) q4(2)],'color',rgb(7,:));
    plot([q4(1) q3(1)],[q4(2) q3(2)],'color',rgb(7,:));
    plot([q3(1) q1(1)],[q3(2) q1(2)],'color',rgb(7,:));
    %--------------------------------------
    plot([q1x(n) x(n)],[q1z(n) z(n)],'k--');
    plot([q2x(n) x(n)],[q2z(n) z(n)],'k--');
    %--------------------------------------
    floorx = [-10 10];
    floorz = [0 0];
    line(floorx,floorz,'Color','green','LineStyle','--');
    
    xlabel('x [m]','Interpreter','Latex'); 
    ylabel('z [m]','Interpreter','Latex'); 
     
    
    axis ([-10 10 0 10]); 
    pbaspect([2 1 1]);
    grid on;
    grid minor;
    title('2D animation','Interpreter','Latex');
    set(findall(gcf,'type','line'),'linewidth',2);
    drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
hold off;
end