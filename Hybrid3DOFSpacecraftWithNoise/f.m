function xdot = f(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Function  
%
% Author: Ricardo Sanfelice and Bharani Malladi
% Revised by: Giulia Zucchini
%
% Project: Simulation of a hybrid system
%
% Name: f.m
%
% Description: Flow map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global xd xd2 Ts 
global K11 K12 K2 K31 K32 K41 K42
global tnoise NOISE2p NOISE3p

xsys = x(1:4);
xsys2 = x(5:6);
h = x(7);
p = x(8);
q = x(9);
tau = x(10);

mu = 3.98600444*10^14; ro = 7100*1000;
n = sqrt(mu/ro^3);
A = [0     0    1 0;
     0     0    0 1;
     3*n^2 0    0 2*n;
     0     0 -2*n 0];
m = 1*500;  
B = [0  0; 0  0; 1/m  0;0  1/m];

A2 = [0     1;
     -n^2  0];
B2 = [0;  1/m];

% new B matrices for phase 4
m4 = 2500;

B4 = [0  0; 0  0; 1/m4  0;0  1/m4];
B42 = [0;  1/m4];

%-------------- CONTROLLER Phase-I --------------% 
if p == 1
    
    %------------ Noise ----------%
    tn = tau;
    tmp = abs(tnoise-tn);
    [row col] = min(tmp);
    rnoise1 = NOISE2p(col)*0.05;
    rnoisev = NOISE2p(col)*0.0001;
    
    %----- position error -----%
    xi(1) = xsys(1) + rnoise1;
    xi(2) = xsys(2) + rnoise1;
    xi(3) = xsys(3) + 0.05*rnoisev;
    xi(4) = xsys(4) + 0.05*rnoisev;
    xi = [xi(1) xi(2) xi(3) xi(4)]';
    
    %------- LQR cotroller for xy system --------%
    up = - K11*(xi-[0; 0; 0; 0]);
    inp = B*up;

    % saturation xy system
    if norm(inp) > 0.02*2/3
        inp = 0.02*2/3* (inp/norm(inp));
    end

    % system for xy component
    xsysdot= A*xsys + inp;
 
    
    %---------- position error -----------%
    xi2(2) = xsys2(1) + 0.001*rnoise1;
    xi2(2) = xsys2(2) + 0.001*0.05*rnoisev;    
    xi2 = [xi2(1) xi2(2)]';
    
    %-------LQR cotroller for z system--------%
    up2 = - K12*(xi2-[0; 0]);
    inp2 = B2*up2;

    % saturation z system
    if norm(inp2) > 0.02*1/3
        inp2 = 0.02*1/3* (inp2/norm(inp2));
    end

    % system for z component
    xsysdot2= A2*xsys2 + inp2;


%-------------- CONTROLLER Phase-II -------------%
elseif p == 2
    
    %------------Noise----------%
    tn = tau;
    tmp = abs(tnoise-tn);
    [row col] = min(tmp);
    rnoise1 = NOISE2p(col)*0.05; % 5% of measurement noise
    rnoisev = NOISE2p(col)*0.0001;
    
    % position error               
    xi(1) = xsys(1) + rnoise1;
    xi(2) = xsys(2) + rnoise1;
    xi(3) = xsys(3) + 0.05*rnoisev;
    xi(4) = xsys(4) + 0.05*rnoisev;

    % angle error
    Th = (atan2(xsys(2),xsys(1)));
    r = sqrt(xi(1)^2+xi(2)^2);
    w =n;
    an = 179;
    % desired angle
    Ths = h*an*2*pi/360;
    k1 = 40;    k2 = 0.1;     k3 = 25;    k4 = 0.047;   
    % Theta original
    vth = (-xi(3)*sin(Th)+xi(4)*cos(Th)); 
    Thd = vth/r;
    
    cosThe = cos(Th)*cos(Ths) - sin(Th)*sin(Ths);
    sinThe = sin(Th)*cos(Ths) - cos(Th)*sin(Ths);

    Te = atan2(sinThe,cosThe);

    vr = (xi(3)*cos(Th)+xi(4)*sin(Th));  %Rho Original
    rd = vr;
    ur = -k1*(rd-0)-k2*(r-150);
    wr = -((3*w^2*xi(1))+xi(4)*(2*w+Thd))*cos(Th)+xi(3)*(2*w+Thd)*sin(Th);  % omegar original
    ar = ur+wr;
    ut = -r*(k3*(Thd-0)+k4*(Te));
    wt = ((3*w^2*xi(1))+xi(4)*(2*w+Thd))*sin(Th)+xi(3)*(2*w+Thd)*cos(Th)+vr*Thd;   %omegatheta original
    at = ut+wt;
    inp = B*[cos(Th) -sin(Th);sin(Th) cos(Th)]*m*[ar;at];    % original input

    % saturation for xy model
    if norm(inp) > 0.02*2/3
        inp = 0.02*2/3* (inp/norm(inp));
    end

    xsysdot = A*xsys + inp;
    
    % position error
    xi2(2) = xsys2(1) + 0.001*rnoise1;
    xi2(2) = xsys2(2) + 0.001*0.05*rnoisev;    
    xi2 = [xi2(1) xi2(2)]';
    
    %-------z cotroller--------%
    up2 = - K2*(xi2-[0; 0]);
    inp2 = B2*up2;

    %saturation for z
    if norm(inp2) > 0.02*1/3
        inp2 = 0.02*1/3* (inp2/norm(inp2));
    end

    %system for z component
    xsysdot2= A2*xsys2 + inp2;

%-------------- CONTROLLER Phase-III -------------%
elseif p == 3
    
    %------------Noise----------%
    tn = tau;
    tmp = abs(tnoise-tn);
    [row col] = min(tmp);
    rnoise2 = NOISE3p(col)*0.05;
    rnoisev = NOISE2p(col)*0.0001;
    
    % position error
    xi(1) = xsys(1)+rnoise2;
    xi(2) = xsys(2)+rnoise2;
    xi(3) = xsys(3)+0.05*rnoisev;
    xi(4) = xsys(4)+0.05*rnoisev;
    xi = [xi(1) xi(2) xi(3) xi(4)]';    
    
    if q == 1
    
        %--- xy LQR controller ---%
        up = -B*K31*(xi-[xd;0;0]);
        inp = up;

    elseif q == 2

        %--- xy P controller ---%
        k_1 = -0.0007;  k_2 = -0.15;    k_3 = -0.006;   k_4 = -0.22;

        ax = 3*n^2*xi(1)+k_1*xi(1) +k_2*(xi(3)-0);
        ay = k_3*xi(2) +k_4*xi(4);
        inp = [0 0 ax ay]';
    
    end 
    
    xsysdot = A*xsys + inp;
    
    % position error
    xi2(2) = xsys2(1) + 0.0001*rnoise2;
    xi2(2) = xsys2(2) + 0.0001*0.05*rnoisev;    
    xi2 = [xi2(1) xi2(2)]';
    
    %--- z LQR controller ---%
    up2 = -B2*K32*(xi2-[xd2;0]);
    inp2 = up2;
    
    xsysdot2 = A2*xsys2 + inp2;   
    
%-------------- CONTROLLER Phase-IV -------------%
elseif p == 4    
    
    % reaching x y position
    xp = [0; 20000];
    % reaching z position
    xp2 = 0;

    %------------Noise----------%
    tn = tau;
    tmp = abs(tnoise-tn);
    [row col] = min(tmp);
    rnoise4 = NOISE2p(col)*0.05;%*0.1;
    rnoisev = NOISE2p(col)*0.0001;%*0.1;

    % position error
    xi(1) = xsys(1)+rnoise4;
    xi(2) = xsys(2)+rnoise4;
    xi(3) = xsys(3)+0.05*rnoisev;
    xi(4) = xsys(4)+0.05*rnoisev;
    xi = [xi(1) xi(2) xi(3) xi(4)]';
    
    %--- xy LQR controller ---%

    up = -B4*K41*(xsys-[xp;0;0]);
    inp = up;

    xsysdot = A*xsys + inp;
    
    % position error
    xi2(2) = xsys2(1) + 0.001*rnoise4;
    xi2(2) = xsys2(2) + 0.001*0.05*rnoisev;    
    xi2 = [xi2(1) xi2(2)]';
    
    %--- z LQR controller ---%
    up2 = -B42*K42*(xi2-[xp2;0]);
    inp2 = up2;

    xsysdot2 = A2*xsys2 + inp2;
    
    
end 
    
hdot = 0;
pdot = 0;
qdot = 0;
taudot = 1;

xdot = [xsysdot;xsysdot2;hdot;pdot;qdot;taudot];

end 
    