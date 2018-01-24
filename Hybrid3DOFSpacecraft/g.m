function xplus = g(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Function  
%
% Author: Ricardo Sanfelice and Bharani Malladi
% Revised by: Giulia Zucchini
%
% Project: Simulation of a hybrid system
%
% Name: g.m
%
% Description: Jump map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global er xd xd2

xsys = x(1:4);
xsys2 = x(5:6);
h = x(7);
p = x(8);
q = x(9);
tau = x(10);

% -------- general switching between phases --------- %
pp = p;
r3 = sqrt(x(1)^2+x(2)^2+x(5)^2);

if (abs(r3-700) <= 0.1) && (p == 1)
    pp = 2;
elseif (abs(r3-150) <= 0.1) && (p == 2)
    pp = 3;
elseif (r3 <= 1e-2) && (p == 3)
    pp = 4;
end 

%--------- phase II ------------ %
hp = h;
Th = (atan2(xsys(2),xsys(1)));
rho = 10*pi/180;

if (Th >= rho && h ==-1) || (Th <= -rho && h ==1)
    hp = -h;
end

%--------- phase III ------------ %
qq = q;

rPhase3 = norm([xsys(1)-xd(1) xsys(2)-xd(2) xsys2(1)-xd2]);

if p == 3 
    if rPhase3<=er && q==1
        qq = 3-q;      
    end
end


xplus =[xsys;xsys2;hp;pp;qq;tau];

end

