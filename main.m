%% Lander with AMEID and control - Numerical caculation
%Date:2021/01/11
%Code writter: Lung-Hui Wu <b07901019@ntu.edu.tw> 
%Professor: Cheng-Wei Chen <cwchenee@ntu.edu.tw>
%All rights reserved.
%----------------------------------------------------------------------------------

clear;clc;
%% initial states setting
qi =[0 11 0.8 0.8 1.05];        %I.C. of generalized coordinates: qi = [z(to) x(to) theta(to) l1(to) l2(to)];
ui =[0 0 0 0 0];          %I.C. of generalized velocities:  ui = [x_dot(to) z_dot(to) l1_dot(to)];
xmi  =[0 0 0];        %I.C. of AMEID states: xmi = [
lenq=length(qi);
%% PID parameters setting-----The digital P-controller is implemented using the so-called velocity form.
An1=10;
An2=10;
Kp1=800;
Kp2=800;
C1 = 0.65;
C2 = 0.65;
R1 = 0.6;
R2 = 0.6;
%% Parameters setting
%% 2-D Lander parameters
mb = 10;       %(kg)- Mass of the main body 
m1 = 0.5;      %(kg)- Mass of footpad 1
m2 = 0.5;      %(kg)- Mass of footpad 2
Ig = 3.7309;   %(kg*m^2) - Moment of inertia of lander
k  = 7000;     %(N/m) - Stiffness of leg
c  = 350;      %(N*s/m) - Damping coefficient of leg
lss = 0.8;    %(m) - Length of passive legs
W = 1;         %(m) - Width of main body
H = 1;         %(m) - Height of main body
R = sqrt(2);   %(m) - Diagonal length of main body
gamma = pi/4;  %(rad) - Gamma
alpha = pi/3;  %(rad) - alpha
S = 1.3856;    %(m) - Arm length to passive leg
%% Gravity & Ground Properties
g = 9.80665;           %(m/s^2)-Standard gravity-the nominal gravitational acceleration of an object in a vacuum near the surface of the Earth
kf = 7.5e4;            %(N/m)-Ground spring stiffness
cf = 130;              %(N.s/m)-Ground damper coefficient
uk = 0.3;              %kinetic friction coefficient
%% AMEID parameter values
md=0.5;                %(kg)
mp=0.2;                %(kg)
La=6.4e-3;             %(H)
Ra=5.2;                %(Ohm)
kv=17;                 %(V.s/m)
kF=17;                 %(N/A)
L_stroke=0.5;          %(m)
%% Ameid state-space matrix
%pair(A,B,C)
Aam = [-Ra/La 0 -kv/La;0 0 1;kF/(mp+md) 0 0];
Bam = [1/La;0;0];
Cam = [kF 0 0];
%% Simulation condition
T = 5;                     %Simulation time 
dt = 0.0005;               %Sampling time
N = floor(T/dt);           %Steps
t = (0:1:N)'*dt;
%% Definition of lander states
qn = zeros(N+1,lenq);      %Generalized coordinates
un = zeros(N+1,lenq);      %Generalized velocities

q1n = zeros(N+1,2);        %Footpad Mass1 coordinates
u1n = zeros(N+1,2);        %Footpad Mass1 velocities

q2n = zeros(N+1,2);        %Footpad Mass1 coordinates
u2n = zeros(N+1,2);        %Footpad Mass1 velocities

F_impactn1 = zeros(N+1,1);  %Footpad impact force1
F_impactn2 = zeros(N+1,1);  %Footpad impact force2
%% Definition of AMEID states
xm1 = zeros(N+1,3);
xm2 = zeros(N+1,3);
%Input
V1 = zeros(N+1,1);
V2 = zeros(N+1,1);
%PID error matrix
ek1 = zeros(N+2,1);
ek2 = zeros(N+2,1);
%AMEID Force matrix
Fmeid1 = zeros(N+1,1);
Fmeid2 = zeros(N+1,1);
%setting initial condition
qn(1,:) = qi;              
un(1,:) = ui;              
xm1(1,:) = xmi;
xm2(1,:) = xmi;
Launch1 = false;
Launch2 = false;

num1 = 0;
count1 = 0;
num2 = 0;
count2 = 0;
tic;

for i = 1:1:N
xbf = qn(i,1);
zbf = qn(i,2);
l1f = qn(i,3);
l2f = qn(i,4);
thetaf = qn(i,5);

xb_dotf = un(i,1);
zb_dotf = un(i,2);
l1_dotf = un(i,3);
l2_dotf = un(i,4);
theta_dotf = un(i,5);

q1n(i,:) = [ xbf - cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) - l1f*cos(alpha)*cos(thetaf) + l1f*sin(alpha)*sin(thetaf),...
             zbf + cos(thetaf)*(H/2 - S*sin(alpha)) - sin(thetaf)*(W/2 + S*cos(alpha)) - l1f*cos(alpha)*sin(thetaf) - l1f*sin(alpha)*cos(thetaf)];

u1n(i,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) - l1_dotf*cos(thetaf)*cos(alpha) + l1_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l1f + theta_dotf*sin(thetaf)*cos(alpha)*l1f ,...
             zb_dotf - theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l1_dotf*cos(thetaf)*sin(alpha) - l1_dotf*sin(thetaf)*cos(alpha) - theta_dotf*cos(thetaf)*cos(alpha)*l1f + theta_dotf*sin(thetaf)*sin(alpha)*l1f];
    
q2n(i,:) = [ xbf + cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) + l2f*cos(alpha)*cos(thetaf) + l2f*sin(alpha)*sin(thetaf),...
             zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)];

u2n(i,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + l2_dotf*cos(thetaf)*cos(alpha) + l2_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l2f - theta_dotf*sin(thetaf)*cos(alpha)*l2f ,...
             zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f];
         
gNb = qn(i,2);
gN1 = q1n(i,2);
gN2 = q2n(i,2);

rNb = un(i,2);
rN1 = u1n(i,2);
rN2 = u2n(i,2);

rTb = un(i,1);
rT1 = u1n(i,1);
rT2 = u2n(i,1);

% Mf in symb
Mn = [m1 + m2 + mb,...
    0,...
    -m1*cos(alpha + thetaf),...
    m2*cos(alpha - thetaf),...
    (m1*(2*l1f*sin(alpha + thetaf) - cos(thetaf)*(H - 2*S*sin(alpha)) + sin(thetaf)*(W + 2*S*cos(alpha))))/2 - (m2*(cos(thetaf)*(H - 2*S*sin(alpha)) + sin(thetaf)*(W + 2*S*cos(alpha)) - 2*l2f*sin(alpha - thetaf)))/2;
    0,...
    m1 + m2 + mb,...
    -m1*sin(alpha + thetaf),...
    -m2*sin(alpha - thetaf),...
    (m2*(cos(thetaf)*(W + 2*S*cos(alpha)) - sin(thetaf)*(H - 2*S*sin(alpha)) + 2*l2f*cos(alpha - thetaf)))/2 - (m1*(2*l1f*cos(alpha + thetaf) + cos(thetaf)*(W + 2*S*cos(alpha)) + sin(thetaf)*(H - 2*S*sin(alpha))))/2;
    -m1*cos(alpha + thetaf),...
    -m1*sin(alpha + thetaf),...
    m1,...
    0,...
    (m1*(H*cos(alpha) + W*sin(alpha)))/2;
    m2*cos(alpha - thetaf),...
    -m2*sin(alpha - thetaf),...
    0,...
    m2,...
    -(m2*(H*cos(alpha) + W*sin(alpha)))/2;
    m1*(l1f*sin(alpha + thetaf) - cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha))) - m2*(cos(thetaf)*(H/2 - S*sin(alpha)) - l2f*sin(alpha - thetaf) + sin(thetaf)*(W/2 + S*cos(alpha))),...
    m2*(l2f*cos(alpha - thetaf) + cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha))) - m1*(l1f*cos(alpha + thetaf) + cos(thetaf)*(W/2 + S*cos(alpha)) + sin(thetaf)*(H/2 - S*sin(alpha))),...
    (m1*(H*cos(alpha) + W*sin(alpha)))/2,...
    -(m2*(H*cos(alpha) + W*sin(alpha)))/2,...
    Ig + (H^2*m1)/4 + (H^2*m2)/4 + S^2*m1 + S^2*m2 + (W^2*m1)/4 + (W^2*m2)/4 + l1f^2*m1 + l2f^2*m2 + 2*S*l1f*m1 + 2*S*l2f*m2 + S*W*m1*cos(alpha) + S*W*m2*cos(alpha) - H*S*m1*sin(alpha) - H*S*m2*sin(alpha) + W*l1f*m1*cos(alpha) + W*l2f*m2*cos(alpha) - H*l1f*m1*sin(alpha) - H*l2f*m2*sin(alpha)
    ];

% ht in symb
hn = [                                       0;
                             -g*(m1 + m2 + mb);
                             g*m1*(cos(alpha)*sin(thetaf) + sin(alpha)*cos(thetaf)) - (k*(2*l1f - 2*lss))/2 - c*l1_dotf;
                             - c*l2_dotf - (k*(2*l2f - 2*lss))/2 - g*m2*(cos(alpha)*sin(thetaf) - sin(alpha)*cos(thetaf));
                             g*m1*(cos(thetaf)*(W/2 + S*cos(alpha)) + sin(thetaf)*(H/2 - S*sin(alpha)) + l1f*cos(alpha)*cos(thetaf) - l1f*sin(alpha)*sin(thetaf)) - g*m2*(cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) + l2f*cos(alpha)*cos(thetaf) + l2f*sin(alpha)*sin(thetaf))];
 
  
LambdaN1 = cf*min(-sign(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf)), 0)*min(-sign(theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - zb_dotf + theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) + l1_dotf*cos(thetaf)*sin(alpha) + l1_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l1f - theta_dotf*sin(thetaf)*sin(alpha)*l1f), 0)*(theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - zb_dotf + theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) + l1_dotf*cos(thetaf)*sin(alpha) + l1_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l1f - theta_dotf*sin(thetaf)*sin(alpha)*l1f) - kf*min(-sign(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf)), 0)*(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf));  
LambdaN2 = kf*min(sign(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)), 0)*(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)) - cf*min(sign(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)), 0)*min(sign(zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f), 0)*(zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f);
LambdaNb = kf*zbf*min(sign(zbf), 0) - cf*zb_dotf*min(sign(zb_dotf), 0)*min(sign(zbf), 0);  
LambdaN = [LambdaN1;LambdaN2;LambdaNb]; 
% 
LambdaT1 = uk*sign(xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) - l1_dotf*cos(thetaf)*cos(alpha) + l1_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l1f + theta_dotf*sin(thetaf)*cos(alpha)*l1f)*(kf*min(-sign(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf)), 0)*(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf)) - cf*min(-sign(sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - zbf + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf)), 0)*min(-sign(theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - zb_dotf + theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) + l1_dotf*cos(thetaf)*sin(alpha) + l1_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l1f - theta_dotf*sin(thetaf)*sin(alpha)*l1f), 0)*(theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - zb_dotf + theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) + l1_dotf*cos(thetaf)*sin(alpha) + l1_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l1f - theta_dotf*sin(thetaf)*sin(alpha)*l1f));
LambdaT2 = -uk*sign(xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + l2_dotf*cos(thetaf)*cos(alpha) + l2_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l2f - theta_dotf*sin(thetaf)*cos(alpha)*l2f)*(kf*min(sign(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)), 0)*(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)) - cf*min(sign(zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)), 0)*min(sign(zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f), 0)*(zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f));
LambdaTb = -uk*sign(xb_dotf)*(kf*zbf*min(sign(zbf), 0) - cf*zb_dotf*min(sign(zb_dotf), 0)*min(sign(zbf), 0));
LambdaT = [LambdaT1;LambdaT2;LambdaTb];
%
WN=[0,0,0;
    1,1,1;
    -cos(alpha)*sin(thetaf)-sin(alpha)*cos(thetaf),0,0;
    0,cos(alpha)*sin(thetaf) - sin(alpha)*cos(thetaf),0;
    l1f*sin(alpha)*sin(thetaf) - sin(thetaf)*(H/2 - S*sin(alpha)) - l1f*cos(alpha)*cos(thetaf) - cos(thetaf)*(W/2 + S*cos(alpha)),cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) + l2f*cos(alpha)*cos(thetaf) + l2f*sin(alpha)*sin(thetaf),0];
WT=[1,1,1;
    0,0,0;
    sin(thetaf)*sin(alpha) - cos(thetaf)*cos(alpha),0,0;
    0,cos(thetaf)*cos(alpha) + sin(thetaf)*sin(alpha),0;
    sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) + cos(thetaf)*sin(alpha)*l1f + sin(thetaf)*cos(alpha)*l1f,cos(thetaf)*sin(alpha)*l2f - sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) - sin(thetaf)*cos(alpha)*l2f,0];
P = WN*LambdaN+WT*LambdaT;
%%%% ------------------AMEID force------------------------%%%%
% left
if Launch1 == false
    ek1(i+2) = qn(i,5);
   V1(i+1) = An1*ek1(i+1)-Kp1*(ek1(i+2)-ek1(i+1));
   if (qn(i,5)>pi*C1 || qn(i,5)<-pi*C1)
       V1(i+1) = V1(i+1) - R1 * An1*ek1(i+1);
   end
   if xm1(i,2)<0
      xm1(i,2)=0;
   end
   xm1(i+1,:) = ((Aam*dt+eye(3))*(xm1(i,:))'+Bam*dt*V1(i+1))';   
   Fl = -Cam*(xm1(i,:))';
   
elseif Launch1 == true
   ek1(i+2) = -(0-u1n(i,2));
   %V1(i+1) = V1(i) + Kp1*(ek1(i+2)-ek1(i+1));
   V1(i+1) = 0;
   
   if xm1(i,2)<0
      xm1(i,2)=0;
   end
   xm1(i+1,:) = ((Aam*dt+eye(3))*(xm1(i,:))'+Bam*dt*V1(i+1))';   
   Fl = -Cam*(xm1(i,:))';
end
% right
if Launch2 == false
    ek2(i+2) = qn(i,5);
   V2(i+1) = -An2*ek2(i+1)-Kp2*(ek2(i+2)-ek2(i+1));
   if (qn(i,5)>pi*C2 || qn(i,5)<-pi*C2)
       V2(i+1) = V2(i+1) + R2 * An2*ek2(i+1);
   end
   if xm2(i,2)<0
      xm2(i,2)=0;
   end
   xm2(i+1,:) = ((Aam*dt+eye(3))*(xm2(i,:))'+Bam*dt*V2(i+1))'; 
   Fr = -Cam*(xm2(i,:))';
elseif Launch2 == true
   ek2(i+2) = -(0-u2n(i,2));
   %V2(i+1) = V2(i) + Kp2*(ek2(i+2)-ek2(i+1));
   V2(i+1) = 0;
   
   if xm2(i,2)<0
      xm2(i,2)=0;
   end
   xm2(i+1,:) = ((Aam*dt+eye(3))*(xm2(i,:))'+Bam*dt*V2(i+1))'; 
   Fr = -Cam*(xm2(i,:))';
end
%-------L_stroke constraint------
if (xm1(i,2) >= L_stroke) || (num1 >= 1)
   xm1(i,2) = L_stroke;
   xm1(i+1,3) = 0;
   Launch1 = false;
   num1 = num1 + 1;
end
if (xm2(i,2) >= L_stroke) || (num2 >= 1)
   xm2(i,2) = L_stroke;
   xm2(i+1,3) = 0;
   Launch = false;
   num2 = num2 + 1;
end
%--------------------------------
hl = [sin(thetaf); -cos(thetaf); 0;0;W+(S+l2f)*cos(alpha)]*Fl;
hr = [sin(thetaf); -cos(thetaf); 0;0;-W-(S+l1f)*cos(alpha)]*Fr;

%---Euler's method (Solve DAE)---
un(i+1,:) = un(i,:) + (Mn\(hn+P+hl+hr)*dt)';
qn(i+1,:) = qn(i,:) + un(i,:)*dt;
%--------------------------------

xbf = qn(i+1,1);
zbf = qn(i+1,2);
l1f = qn(i+1,3);
l2f = qn(i+1,4);
thetaf = qn(i+1,5);

xb_dotf = un(i+1,1);
zb_dotf = un(i+1,2);
l1_dotf = un(i+1,3);
l2_dotf = un(i+1,4);
theta_dotf = un(i+1,5);

u1n(i+1,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) - l1_dotf*cos(thetaf)*cos(alpha) + l1_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l1f + theta_dotf*sin(thetaf)*cos(alpha)*l1f ,...
             zb_dotf - theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l1_dotf*cos(thetaf)*sin(alpha) - l1_dotf*sin(thetaf)*cos(alpha) - theta_dotf*cos(thetaf)*cos(alpha)*l1f + theta_dotf*sin(thetaf)*sin(alpha)*l1f];

u2n(i+1,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + l2_dotf*cos(thetaf)*cos(alpha) + l2_dotf*sin(thetaf)*sin(alpha) + theta_dotf*cos(thetaf)*sin(alpha)*l2f - theta_dotf*sin(thetaf)*cos(alpha)*l2f ,...
             zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l2_dotf*cos(thetaf)*sin(alpha) + l2_dotf*sin(thetaf)*cos(alpha) + theta_dotf*cos(thetaf)*cos(alpha)*l2f + theta_dotf*sin(thetaf)*sin(alpha)*l2f];
 
                    
%-----Does Lander contact to the ground?---------           
if (u1n(i+1,2)*u1n(i,2)<0)
   if (xm1(i,2)<L_stroke)
      Launch1 = true;
   else
      Launch1 = false;
   end
end
if (u2n(i+1,2)*u2n(i,2)<0)
   if (xm2(i,2)<L_stroke)
      Launch2 = true;
   else
      Launch2 = false;
   end
end
if (qn(i,2)>2.5)
    Launch1 = false;
    Launch2 = false;
end
F_impactn1(i) = -k*(qn(i,3)-lss)-c*un(i,3);
F_impactn2(i) = -k*(qn(i,4)-lss)-c*un(i,4);

end
xm1(end,2) = xm1(end-1,2);
xm2(end,2) = xm2(end-1,2);
Fmax=max(F_impactn1,F_impactn2);

toc;

if abs(un(end,2))<=1e-2
    flag='stable';
else
    flag='unstable';
end





%% Plot non-linear result
figure(1);
subplot(2,1,1);
plot(t,qn(:,2)); %zbf
L0=legend('$z_{b}$f [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Height(m)','Interpreter','Latex');
title({'1D Lander Free Fall Landing';['xb(t0)=',num2str(qn(1,1)),'[m],','zb(t0)=',num2str(qn(1,2)),'[m],',...
    '$l_{1}$(t0)=',num2str(qn(1,3)),'[m]','$l_{2}$(t0)=',num2str(qn(1,4)),'[m]','theta(t0)=',num2str(qn(1,5)),'[rad]'];...
    ['  $\dot{xb}$(t0)=',num2str(un(1,1)),'[m/s],','$\dot{zb}$(t0)=',num2str(un(1,2)),'[m/s],','$\dot{l_{1}}$(t0)=',num2str(un(1,3)),'[m/s].','$\dot{l_{2}}$(t0)=',num2str(un(1,4)),'[m/s].','$\dot{theta}$(t0)=',num2str(un(1,5)),'[rad/s].']},'Interpreter','Latex');
%set(findall(gcf,'type','line'),'linewidth',2);
%kfig(2)=figure;
subplot(2,1,2);
plot(t,qn(:,5));
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Angle(rad)','Interpreter','Latex');
title('Theta v.s. time','Interpreter','Latex');
sgtitle(['\bf Lander Landing status : ',flag],'FontSize',10,'Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);

%{
figure(1);
subplot(3,1,1);
plot(t,qn(:,2)); %zbf
hold on;
%plot(t,qn(:,1)); %xbf
plot(t,qn(:,3)); %l1f
plot(t,qn(:,4)); %l2f
plot(t,qn(:,5)); %theta
hold off;
L0=legend('$z_{b}$f [m]','$l_{1}$f [m]','$l_{2}$f [m]','thetaf [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Response','Interpreter','Latex');
title({'1D Lander Free Fall Landing';['xb(t0)=',num2str(qn(1,1)),'[m],','zb(t0)=',num2str(qn(1,2)),'[m],',...
    '$l_{1}$(t0)=',num2str(qn(1,3)),'[m]'];['  $\dot{xb}$(t0)=',num2str(un(1,1)),'[m/s],','$\dot{zb}$(t0)=',num2str(un(1,2)),'[m/s],','$\dot{l_{1}}$(t0)=',num2str(un(1,3)),'[m/s].']},'Interpreter','Latex');
%set(findall(gcf,'type','line'),'linewidth',2);
%kfig(2)=figure;
subplot(3,1,2);
plot(t,F_impactn1);hold on;
plot(t,F_impactn2);hold off;
L0=legend('F1 [N](Impact Force on The Footpad)','F2 [N](Impact Force on The Footpad)','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Impact Force','Interpreter','Latex');
title(['Fmax:','idk','[N]'],'Interpreter','Latex');
subplot(3,1,3);
plot(t,un(:,2)); %zb_dotf
L0=legend('$\dot{zb}$f [m/s]','Interpreter','Latex','Location','Southeast');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Velocity Response','Interpreter','Latex');
sgtitle(['\bf Lander Landing status : ',flag],'FontSize',10,'Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);
%}

%{
figure(2);
plot(t,V1);hold on;
plot(t,V2);hold off;
L0=legend('$V$f [V]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;

xlabel('Time [s]','Interpreter','Latex');
ylabel('Voltage','Interpreter','Latex');
title('Control Effort','Interpreter','Latex');
set(findall(gcf,'type','stair'),'linewidth',2);

figure(3);
plot(t,xm1(:,2));hold on;
plot(t,xm2(:,2));hold off;
L0=legend('$xp$f [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;

xlabel('Time [s]','Interpreter','Latex');
ylabel('VCM Plate Position [m]','Interpreter','Latex');
title('AMEID-VCM Plate Position','Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);
%}
%{
%% -------------Animation----------------
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
figure(4);
axis ([-10 10 -0.1 7]); 
axis manual;

for n = 1:40:N
    
    Rc=[cos(qn(n,5)),-sin(qn(n,5));sin(qn(n,5)),cos(qn(n,5))];
    qb=[qn(n,1);qn(n,2)];
    
    r1 = [-W/2;H/2];
    r2 = [W/2;H/2]; 
    r3 = [-W/2;-H/2];
    r4 = [W/2;-H/2];
    
    r1_c = Rc*r1;
    r2_c = Rc*r2;
    r3_c = Rc*r3;
    r4_c = Rc*r4;
    
    q1 = qb + r1_c;
    q2 = qb + r2_c;
    q3 = qb + r3_c;
    q4 = qb + r4_c;
    
    plot(qn(n,1),qn(n,2),'ro');
    hold on;
    plot(q1n(n,1),q1n(n,2),'bs');
    plot(q2n(n,1),q2n(n,2),'bs');
    %--------------Lander body-------------
    plot([q1(1) q2(1)],[q1(2) q2(2)],'color',rgb(7,:));
    plot([q2(1) q4(1)],[q2(2) q4(2)],'color',rgb(7,:));
    plot([q4(1) q3(1)],[q4(2) q3(2)],'color',rgb(7,:));
    plot([q3(1) q1(1)],[q3(2) q1(2)],'color',rgb(7,:));
    %--------------------------------------
    plot([q1n(n,1) q1(1)],[q1n(n,2) q1(2)],'k:');
    plot([q2n(n,1) q2(1)],[q2n(n,2) q2(2)],'k:');
    %--------------------------------------
    hold off;
    floorx = [-10 10];
    floorz = [0 0];
    line(floorx,floorz,'Color','green','LineStyle','--');
    
    xlabel('x [m]','Interpreter','Latex'); 
    ylabel('z [m]','Interpreter','Latex'); 
     
    
    axis ([-10 10 -0.1 11]); 
    pbaspect([200 110 1]);
    grid on;
    grid minor;
    title('2D animation','Interpreter','Latex');
    set(findall(gcf,'type','line'),'linewidth',2);
    drawnow 
    
    if n==1
        pause(2);
    end
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
%}
