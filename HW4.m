%% HW #4
% Submitted by Jesse Austin Stringfellow,    
% Due Oct. 2, 2019
%% Problem #1
%a) See attached SE2.m code
%b) Verification Code
g1 = SE2([1;2], pi/3);
g2 = SE2([-2;1], pi/6);
g3 = adjoint(g1, g2)
%% Problem #2
%From HW #1
gOA = SE2([5;12], pi/3); 
gOB = SE2([2;-1], pi);
%Given:
h = SE2([1,0],pi/2);
%a) Compute hammer head configurations A'(Ap) and B'(Bp) relative to O
gOAp = gOA*h
gOBp = gOB*h
%b) Compoute gA'B' (gApBp).
%gA'B' = gA'O*gOB' where gA'O = inv(gOA')
gApBp = inv(gOAp)*gOBp
%c)
% The idea of this problem is that gApBp can be found by using the Adjoint
% Basically gA'B' = (gA'O)*gOB' and from part a) it can be seen that 
% gOA' = gOA*h and gOB' = gOB*h so
% gA'B' = inv(gOA*h)*(gOB*h) = inv(h) * gAO * gOB * h = inv(h)* gAB * h
% This form can be seen to be the same as the Adjoint(inv(h),gAB)
gAB = inv(gOA)*gOB;
gApBp_2 = adjoint(inv(h),gAB)
%% Problem #3
% Given
gH1H2 = SE2([3.151;3.906],-7*pi/12); % Movement from H1 (Position 1) to H2 (position 2)
gOH2 = SE2([3.243;2.512],[-.131,-.991;.991,-.131]); % Position 2 of end effector relative to base, O.
gHK = SE2([3;-1],-pi/4); %Key relative to end effector frame
%a) question is asking to compute gOH1
%gOH1 = gOH2*gH2H1 where gH2H1 = inv(gH1H2)
gOH1 = gOH2*inv(gH1H2)
%b)compute gOK in positon 2
gOK = gOH2*gHK
%c)Compute gK1K2, which is the movement of the key frame from Position 1 to
%to positon 2
gK1K2 = adjoint(inv(gHK),gH1H2)
%This gives gKH * gH1H2 * gHK and with the leftward order of operations
%'cancels out' the H1

%% Problem #4
%a)
theta = pi/4;
Rx = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)]
%b)
thetay = pi/6; thetaz = 3*pi/4;
Ry = [cos(thetay) 0 sin(thetay);0 1 0; -sin(thetay) 0 cos(thetay)];
Rz = [cos(thetaz) -sin(thetaz) 0;sin(thetaz) cos(thetaz) 0 ;0 0 1];
Ry*Rz
%c)
thetax = pi/3; thetaz = -pi/4;
Rx = [1 0 0;0 cos(thetax) -sin(thetax);0 sin(thetax) cos(thetax)];
Rz = [cos(thetaz) -sin(thetaz) 0;sin(thetaz) cos(thetaz) 0 ;0 0 1];
Rx*Rz
%d) Show that Ry(theta_1)*Ry(theta_2) = Ry(theta_1 + theta_2)







%



