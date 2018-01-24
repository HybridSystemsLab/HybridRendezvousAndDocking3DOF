function v  = C(x) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab Function  
%
% Author: Ricardo Sanfelice and Bharani Malladi
% Revised by: Giulia Zucchini
%
% Project: Simulation of a hybrid system
%
% Name: C.m
%
% Description: Flow set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global er xd xd2

xsys = x(1:4);
xsys2 = x(5:6);
h = x(7);
p = x(8);


v = 0;

% ----- continuing flowing in phase I or phase II ------ %
r3 = sqrt(x(1)^2+x(2)^2+x(5)^2);

if (abs(r3-700) >= 0.1 && p == 1) || (abs(r3-700) <= 0.1 && p == 2) || (abs(r3-150) <= 0.1 &&  p == 3) || (r3 >= 1e-2 && p==4)
    v = 1;
end

% ------ internal condition for flow in phase II ------ %
Th = (atan2(xsys(2),xsys(1)));
rho = 10*pi/180;

if (Th <= rho && h == -1) || (Th >= -rho && h == 1)
    v = 1;
end

% ------ internal flow condition for phase III ----- %
rPhase3 = norm([xsys(1)-xd(1) xsys(2)-xd(2) xsys2(1)-xd2]);

if p == 3 && (rPhase3 >= er)
    v = 1;
end

end

