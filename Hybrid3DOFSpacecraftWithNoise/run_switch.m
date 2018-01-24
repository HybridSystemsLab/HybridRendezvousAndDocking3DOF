format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Main file  
%
% Author: Ricardo Sanfelice and Bharani Malladi
% Revised by: Giulia Zucchini
%
% Project: Simulation of a hybrid system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global er xd xd2 Ts 
global K11 K12 K2 K31 K32 K41 K42
global tnoise NOISE2p NOISE3p 

er = 1;                 % unit distance 
xd = [-25;0]; xd2 = 0;  % desired position to reach in phase III
Ts = 10;

mu = 3.98600444*10^14; ro = 7100*1000;
n = sqrt(mu/ro^3); 
A = [0     0    1   0;
     0     0    0   1;
     3*n^2 0    0   2*n;
     0     0 -2*n   0];
m = 1*500;  
B = [0  0; 0  0; 1/m  0;0  1/m];

A2 = [0     1;
     -n^2  0];
B2 = [0;  1/m];

% new B matrices for phase 4
m4 = 2500;

B4 = [0  0; 0  0; 1/m4  0;0  1/m4];
B42 = [0;  1/m4];

%-----------------NOISE-----------%

N=150000;
Fs = 1000;
tnoise = (0:N-1)/Fs;
sigma = 1;
NOISE2p = (10)*sigma*randn(size(tnoise));
NOISE3p = (1/100)*sigma*randn(size(tnoise));


% ------- Phase I ---------- %

% LQR xy system
Q11 = .1*(1.5e-1)*eye(4);
R11 = [20e4 0; 0 11e4];
[K11,s,e] = lqr(A,B,Q11,R11); 

% LQR for z controller 
Q12 = .1*(1.5e-1)*eye(2);
R12 = 99e3;
 [K12,s,e] = lqr(A2,B2,Q12,R12); 


% --------- Phase II --------%

% LQR for z controller 
Q2 = [138 0;0 10];        R2 = 30e6;
[K2,s,e] = lqr(A2,B2,Q2,R2);


% --------- Phase III --------%

% LQR xy system q = 1 
Q31 = 38.4*eye(4);
R31 = 9.7e3*eye(2);
[K31,s,e] = lqr(A,B,Q31,R31);

% LQR z system q = 1 
Q32 = [138 0;0 10];       
R32 = 30e6;
[K32,s,e] = lqr(A2,B2,Q32,R32);

% --------- Phase IV ----------- %

% LQR for xy system
Q41 = 6e-1*eye(4);
R41 = 11e4*eye(2);
 [K41,s,e] = lqr(A,B4,Q41,R41);

% LQR for z system 
% Q42 = eye(2);
% R42 = 1;
Q42 = [138 0;0 10];       
R42 = 30e6;
[K42,s,e] = lqr(A2,B42,Q42,R42);

% --------- Initial conditions ---------- %
% xi(:,1:2) position in xy plane; xi(:,3:4) velocity in xy plane
% zi(:,1:2) position and velocity in z plane

    % starting from phase 1
 xi = [-2453.7  0  0.5  0];  zi = [6550.19  0.5];           
% xi = [0  6550.19  0  0.5];    zi = [-2453.7  0.5];         
% xi = [8660.3 5000 .5 .5];   zi = [0 0];                    
% xi = [8660.3  0  0.5  0];    zi = [5000  0.5];            

    % starting from phase I far away
% xi = [0 -10000 0.5 0.5];    zi = [0 0];                   
% xi = [-7071 7071 .5 .5];    zi = [0 0];                   
% xi = [-2453.7  0  0.5  0];  zi = [6550.19  0.5];           
% xi = [5000  0 .5  0 ];  zi = [8660.3  .5];                 
% xi = [0  0 0 0.5];  zi = [10000  0.5];                     
% xi = [0  -2453.7  0   0.5 ];  zi = [ 6550.19  0.5];       
% xi = [0 5000  0  .5];  zi = [8660.3  .5];                 

% h is the logic state variable of the controller in Phase II, can be 1 or -1
hint = 1;  %or-1
% define the phase 
pint = 1;  
% q is required in phase III to switch between controllers
qint = 1;
% time variable
tauint = 0;

xint = [xi zi hint pint qint tauint]; 

% initial condition
x0 = xint;

% simulation horizon
T = [0 32000];                                                                
J = [0 10000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t j x] = HyEQsolver( @f,@g,@C,@D,x0',T,J,rule,options);

%% p
figure(1)
plot(t,x(:,8),'r')
hold on 
grid on
xlabel('t')                    
ylabel('phase')


%% 3D MIMO PLOT
figure(10)
plot3(-x(:,2)/1000,x(:,1)/1000,x(:,5)/1000,'b'); 
hold on
grid on
xlabel('Y')                    
ylabel('X')
zlabel('Z')


%% plot hybrid arc _ position
figure(11)
plot(t,x(:,1)/1000,'b'); 
hold on
plot(t,x(:,2)/1000,'r');
grid on
hold on
plot(t,x(:,5)/1000,'g'); 
hold on
xlabel('time (sec)')                    
ylabel('x,y,z (KM)')  
legend('x','y','z')

%% plot hybrid arc _ velocity
figure(12)
plot(t,x(:,3)/1000,'b'); 
hold on
plot(t,x(:,4)/1000,'r');
grid on
hold on
plot(t,x(:,6)/1000,'g'); 
hold on
xlabel('time (sec)')                    
ylabel('$\dot x$,$\dot y$,$\dot z$ (KM/sec)')                      
legend('$\dot x$','$\dot y$','$\dot z$')

%% 3D plot
figure(20)
% plot sphere
view(73,10)
plot3(0,0,0,'b*','linewidth',2)
hold on
% bigger sphere 10 km
% [x1, y1, z1] = sphere;
%  plot3(10*x1, 10*y1, 10*z1,'color',[0 0 0],'LineStyle',':'); %set the radius 10km
%  hold on
%  contour3(10*x1, 10*y1, 10*z1,'k:');
%  hold on

% intermediate sphere 1 km
[x2, y2, z2] = sphere;
%  plot3(1*x2, 1*y2, 1*z2,'color',[0 0 0],'LineStyle',':'); %set the radius 1 km
%  hold on
%  contour3(1*x2, 1*y2, 1*z2,'g:');
plot3(1*x2, 1*y2, 1*z2,'k:');
 hold on
 
 % smaller sphere 0.1 km
[x3, y3, z3] = sphere;
%  plot3(0.1*x3, 0.1*y3, 0.1*z3,'color',[0 0 0],'LineStyle',':'); %set the radius
%  hold on
%  contour3(0.1*x3, 0.1*y3, 0.1*z3,'m:'); 
plot3(0.1*x3, 0.1*y3, 0.1*z3,'m:');
 hold on
grid on
plot3(x(:,2)/1000,x(:,5)/1000,x(:,1)/1000,'r','linewidth',2)
xlabel('Local horizontal - y axis')
h=get(gca,'xlabel');
set(h,'rotation',-28)
zlabel('Local verical - x axis')
ylabel('Out of plane - z axis')
h=get(gca,'ylabel');
set(h,'rotation',0)
hold on

%plot cone
r = -linspace(0,1,10);
theta = linspace(0,2*pi,30);
[r,theta] = meshgrid(r,theta);
conex = r;
coney = r.*cos(theta);
conez = r.*sin(theta);
plot3(-coney,-conez,conex,'c');
hold on 
contour3(-coney,conez,conex,':');
hold on
plot3(20,0,0,'k o','linewidth',3)
hold on
