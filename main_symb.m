%2-D spacecraft lander modeling
%109-1 Control System
%Code writter: Lung-Hui Wu <b07901019@ntu.edu.tw> 
%Professor: Cheng-Wei Chen <cwchenee@ntu.edu.tw>
%All rights reserved.
%----------------------------------------------------------------------------------

clear;clc;
%% -----------Start of the dynamics-----------
%% kinetic energy function
% function yke=ke(m,u) %kinetic energy
% L=length(u);
% yke=0;
% for i=1:L
%     yke=yke+0.5*m*(u(i))^2;
% end
%% Gravity potential energy function
% function ygpe=gpe(m,g,h) %Gravity potential energy
% ygpe=m*g*h;
%% Spring potential energy function
% function yspe=spe(k,x)   %Spring potential energy
% yspe=0.5*k*(x^2);
%-------------------------------------------------

%parameters
syms mb m1 m2 Ig k c lss
syms W H R gamma alpha S

syms g kf cf uk

%% variables
syms xbf zbf l1f l2f thetaf
syms xb_dotf zb_dotf l1_dotf l2_dotf theta_dotf
%replaced by
syms xb(t) zb(t) l1(t) l2(t) theta(t)
syms xb_dot(t) zb_dot(t) l1_dot(t) l2_dot(t) theta_dot(t)

xb_dot(t) = diff(xb,t);
zb_dot(t) = diff(zb,t);
l1_dot(t) = diff(l1,t);
l2_dot(t) = diff(l2,t);
theta_dot(t) = diff(theta,t);

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


%% Generalized coordinates q=[xb,zb,l1,l2,theta]'
q = [xbf;zbf;l1f;l2f;thetaf];
%% q(t)
qt = subs(q,{xbf,zbf,l1f,l2f,thetaf},{xb,zb,l1,l2,theta});
%% Generalized velocities u=dq/dt=[xb_dot,zb_dot,l1_dot,l2_dot,theta_dot]'
u = [xb_dotf;zb_dotf;l1_dotf;l2_dotf;theta_dotf];
%% u(t)
ut = subs(u,{xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf},{xb_dot,zb_dot,l1_dot,l2_dot,theta_dot});

lenq = length(q);   %don't lnow if it is necessary?
%% vectors
qb = [xbf;zbf]; 
ub = [xb_dotf;zb_dotf];
ra1c = [S*cos(alpha)+W/2;S*sin(alpha)-H/2];
ra2c = [-S*cos(alpha)-W/2;S*sin(alpha)-H/2];
ql1c = [l1f*cos(alpha);l1f*sin(alpha)];
ql2c = [-l2f*cos(alpha);l2f*sin(alpha)];

Rc = [cos(thetaf) -sin(thetaf);sin(thetaf) cos(thetaf)];

ra1 = Rc * ra1c;
ra2 = Rc * ra2c;
ql1 = Rc * ql1c;
ql2 = Rc * ql2c;

qs1f = qb - ra1;
qs2f = qb - ra2;
q1f = qs1f - ql1;
q2f = qs2f - ql2;

qs1 = subs(qs1f,{xbf, zbf, thetaf},{xb, zb, theta});
qs2 = subs(qs2f,{xbf, zbf, thetaf},{xb, zb, theta});
q1 = subs(q1f,{xbf, zbf, thetaf, l1f},{xb, zb, theta, l1});
q2 = subs(q2f,{xbf, zbf, thetaf, l2f},{xb, zb, theta, l2});
u1 = diff(q1, t);
u2 = diff(q2, t);

u1f = subs(u1,{diff(xb(t),t),diff(zb(t),t),diff(l1(t),t),diff(l2(t),t),diff(theta(t),t)},{xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf});
u2f = subs(u2,{diff(xb(t),t),diff(zb(t),t),diff(l1(t),t),diff(l2(t),t),diff(theta(t),t)},{xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf});
% Kinetic Energy
Tbf = ke(mb,ub);
T1f = ke(m1,u1f);
T2f = ke(m2,u2f);
Tig = ke(Ig,theta_dotf);
Tf = Tbf+T1f+T2f+Tig;
% Potential Energy /gpe: gravity potential, spe: spring potential
Vf = gpe(mb,g,qb(2)) + gpe(m1,g,q1f(2)) + gpe(m2,g,q2f(2)) + spe(k,l1f-lss) + spe(k,l2f-lss);
% Dissipative Energy
Wf = spe(c,l1_dotf) + spe(c,l2_dotf);
%Lagrangian
Lf = Tf-Vf;
Qd = [diff(Wf,u(1));diff(Wf,u(2));diff(Wf,u(3));diff(Wf,u(4));diff(Wf,u(5))];
L1 = [diff(Lf,u(1));diff(Lf,u(2));diff(Lf,u(3));diff(Lf,u(4));diff(Lf,u(5))];
L2 = diff(L1,t);
L3 = [diff(Lf,q(1));diff(Lf,q(2));diff(Lf,q(3));diff(Lf,q(4));diff(Lf,q(5))];
Qcf = L2 - L3 + Qd;
Qc = subs(Qcf,{xbf,zbf,l1f,l2f,thetaf,xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf},{xb,zb,l1,l2,theta,xb_dot,zb_dot,l1_dot,l2_dot,theta_dot});
%M
Mf1 = [diff(Tf,u(1));diff(Tf,u(2));diff(Tf,u(3));diff(Tf,u(4));diff(Tf,u(5))];
Mf = [diff(Mf1,u(1)) diff(Mf1,u(2)) diff(Mf1,u(3)) diff(Mf1,u(4)) diff(Mf1,u(5))];
Mf = simplify(Mf);
%M^(-1)
Mf_inv = inv(Mf);
Mf_invt=subs(Mf_inv,{xbf,zbf,l1f,l2f,thetaf,xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf},{xb,zb,l1,l2,theta,xb_dot,zb_dot,l1_dot,l2_dot,theta_dot});
Mf_invt = simplify(Mf_invt);
%% Force matrix
h0 = [diff(Wf,u(1));diff(Wf,u(2));diff(Wf,u(3));diff(Wf,u(4));diff(Wf,u(5))];
h1 = [diff(Tf,q(1));diff(Tf,q(2));diff(Tf,q(3));diff(Tf,q(4));diff(Tf,q(5))];   %(dTf/dq)'
h2 = [diff(Vf,q(1));diff(Vf,q(2));diff(Vf,q(3));diff(Vf,q(4));diff(Vf,q(5))];   %(dVf/dq)'
h3 = [diff(h1,u(1)) diff(h1,u(2)) diff(h1,u(3)) diff(h1,u(4)) diff(h1,u(5))]*u; %(ddTf/dudq)'*u

hf=-h0+h1-h2-h3;
ht=subs(hf,{xbf,zbf,l1f,l2f,thetaf,xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf},{xb,zb,l1,l2,theta,xb_dot,zb_dot,l1_dot,l2_dot,theta_dot});
ht=simplify(ht);
%% accelaeration
u_dot = Mf_invt * (ht + Qc);
%------------------------------------------   
%% minimize symbolic function
symmin=@(x,y)feval(symengine,'min',x,y);
%---------------------------------------------
%% Contact force
%contact distance (normal)
syms gNb gN1 gN2       %zb(t) q1(2)
%contact velocity (normal)
syms rNb rN1 rN2       %zb_dot(t) u1(2)
%contact velocity (tangential)
syms rTb rT1 rT2      %xb_dot(t) u1(1)

syms kf cf         %ground spring and damper
syms uk            %kinetic friction coefficient
gNb = q(2);
gN1 = q1f(2);
gN2 = q2f(2);
gN = [gN1 gN2 gNb];

rNb = u(2);
rN1 = u1f(2);
rN2 = u2f(2);
rN = [rN1 rN2 rNb];

rTb = u(1);
rT1 = u1f(1);
rT2 = u2f(1);
rT = [rT1 rT2 rTb];
%{
%% Gravity & Ground Properties
g = 9.80665;             %(m/s^2)-Standard gravity-the nominal gravitational acceleration of an object in a vacuum near the surface of the Earth
kf = 7.5e4;              %(N/m)-Ground spring stiffness
cf = 130;                %(N.s/m)-Ground damper coefficient
uk = 0.3;                %kinetic friction coefficient
%}

LambdaN1 = kf*gN1*symmin(sign(gN1),0) - cf*rN1*symmin(sign(rN1),0)*symmin(sign(gN1),0);  
LambdaN2 = kf*gN2*symmin(sign(gN2),0) - cf*rN2*symmin(sign(rN2),0)*symmin(sign(gN2),0);  
LambdaNb = kf*gNb*symmin(sign(gNb),0) - cf*rNb*symmin(sign(rNb),0)*symmin(sign(gNb),0);
LambdaN = [LambdaN1;LambdaN2;LambdaNb];
% 
LambdaT1 = -uk*LambdaN1*sign(rT1);
LambdaT2 = -uk*LambdaN2*sign(rT2);
LambdaTb = -uk*LambdaNb*sign(rTb);
LambdaT = [LambdaT1;LambdaT2;LambdaTb];
%
WN = [diff(gN,q(1));diff(gN,q(2));diff(gN,q(3));diff(gN,q(4));diff(gN,q(5))];
WT = [diff(rT,u(1));diff(rT,u(2));diff(rT,u(3));diff(rT,u(4));diff(rT,u(5))];
P = WN*LambdaN+WT*LambdaT;

Pt=subs(P,{xbf,zbf,l1f,l2f,thetaf,xb_dotf,zb_dotf,l1_dotf,l2_dotf,theta_dotf},{xb,zb,l1,l2,theta,xb_dot,zb_dot,l1_dot,l2_dot,theta_dot});

%--------------------------------
%% Without Ameid
M_bar = Mf_invt*(Pt+ht);  %Withot Ameid
%--------------------------------
%% With Ameid 
%AMEID dynamic
%parameters
syms kF mp md La kv Ra kd L_stroke
%states
syms iarf xprf xpr_dotf 
syms iar(t) xpr(t) xpr_dot(t)     %xmeidrt(t)=[ia(t) xp(t) xp_dot(t)]';
syms ialf xplf xpl_dotf 
syms ial(t) xpl(t) xpl_dot(t)     %xmeidlt(t)=[ia(t) xp(t) xp_dot(t)]';
%input
syms Vr Vl dt

xpr_dot(t)=diff(xpr,t);
xpl_dot(t)=diff(xpl,t);
xmeidrf = [iarf;xprf;xpr_dotf];
xmeidrt = subs(xmeidrf,{iarf,xprf,xpr_dotf},{iar,xpr,xpr_dot});
xmeidlf = [ialf;xplf;xpl_dotf];
xmeidlt = subs(xmeidlf,{ialf,xplf,xpl_dotf},{ial,xpl,xpl_dot});

%parameter values
md=0.5;                %(kg)
mp=0.2;                %(kg)
La=6.4e-3;             %(H)
Ra=5.2;                %(Ohm)
kv=17;                 %(V.s/m)
kF=17;                 %(N/A)
kd=0;                  %(N/m)
L_stroke=0.5;          %(m)


%% Ameid state-space matrix
Aam = [-Ra/La 0 -kv/La;0 0 1;kF/(mp+md) 0 0];
Bam = [1/La;0;0];
Cam = [kF 0 0];

Fmeidr = Cam*xmeidrf;
Fmeidl = Cam*xmeidlf;
%% AMEID force matrix
hrff = [sin(thetaf);-cos(thetaf);sym(0);sym(0);-W-(S+l1f)*cos(alpha)];
hrf = hrff * Fmeidr;
hlff = [sin(thetaf);-cos(thetaf);sym(0);sym(0);W+(S+l2f)*cos(alpha)];
hlf = hlff * Fmeidl;
hr = subs(hrf,{xmeidrf(1)},{xmeidrt(1)});
hl = subs(hlf,{xmeidlf(1)},{xmeidlt(1)});
%--------------------------------
M_bar_u=Mf_invt*(Pt+ht+hr+hl); %With AMEID   %dimension: 3x1

%{
%%----------------------------------------------------------------------------------
%%----------------------------------------------------------------------------------
%% Non-linear ode
%non-linear odes
ode1 = diff(ut(1)) == M_bar(1);
ode2 = diff(ut(2)) == M_bar(2);
ode3 = diff(ut(3)) == M_bar(3);
ode4 = diff(ut(4)) == M_bar(4);
ode5 = diff(ut(5)) == M_bar(5);
odes=[ode1;ode2;ode3;ode4;ode5];

%% non-linear fuction
[Veq,Seq] = odeToVectorField(odes);
Meq = matlabFunction(Veq,'vars', {'t','Y'});


%% Simuation calculation
Tf=3;  %simulation time
tspan=[0 Tf]; %simulation horizon
dt=.0005; %time step
x0=[11 0 0 0 0.5 0 0.8 0 0.8 0];  %initial condition
%options = odeset('Mass',eye(10),'RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);
[t,sol] = ode45(Meq,tspan,x0); %use ode45,ode23,ode15s solver to solve the equation
t0 = (0:dt:Tf)';
sol0 = interp1(t,sol,t0,'makima');   %interpolation for equal dt

F_impact= -k.*(sol0(:,1)-lss)-c.*(sol0(:,2));
Fmax=max(F_impact);

%% Plot non-linear result
kfig(1)=figure;
subplot(3,1,1);
plot(t0,sol0(:,7)); % l1(t)
hold on;
plot(t0,sol0(:,9)); % l2(t)
hold off;
L0=legend('$l_{1}$(t)','$l_{2}$(t)','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Response [m]','Interpreter','Latex');
title({'Response without AMEID';['$x_{b}$(t0)=',num2str(x0(3)),'[m],','$z_{b}$(t0)=',num2str(x0(1)),'[m],','$theta$(t0)=',num2str(x0(5)),'[m],',...
    '$l_{1}$(t0)=',num2str(x0(7)),'[m]','$l_{2}$(t0)=',num2str(x0(9)),'[m]'];['  $\dot{xb}$(t0)=',num2str(x0(4)),'[m/s],','$\dot{zb}$(t0)=',num2str(x0(2)),'[m/s],','$\dot{theta}$(t0)=',num2str(x0(6)),'[m/s].','$\dot{l_{1}}$(t0)=',num2str(x0(8)),'[m/s].','$\dot{l_{2}}$(t0)=',num2str(x0(10)),'[m/s].']},'Interpreter','Latex');
%set(findall(gcf,'type','line'),'linewidth',2);
%kfig(2)=figure;
subplot(3,1,2);
plot(t0,sol0(:,5)); % theta(t)
L0=legend('Theta','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('[rad]','Interpreter','Latex');
title('Theta','Interpreter','Latex');
%set(findall(gcf,'type','line'),'linewidth',2);
%kfig(3)=figure;
subplot(3,1,3);
plot(t0,sol0(:,1)); % zb
hold on;
plot(t0,sol0(:,2)); % zb_dot(t)
hold off;
L0=legend('$z_{b}$(t)','$\dot{zb}$(t)','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('zb[m] / zb_dot[m/s]','Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);
%}
%}