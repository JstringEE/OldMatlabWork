%% HW #5
% Submitted by Jesse Austin Stringfellow,    
% Due Oct. 23, 2019
%% Problem #1
%Given
gOA = SE2([4;2], pi/4);
gOB = SE2([5;4], pi/2);
gBC = SE2([-2;0], -pi/6);
% a)
zA1 = [3;6;-pi/12];
% zB1 = gBA*zA1 where gBA = gBO*gOA = inv(gOB)*gOA
gBA = inv(gOB)*gOA;
%gBA.M(1:2,3) = 0;
%The form for g*z  = [R(theta) 0; 0 1] * [z1;z2;z3]. Thus the displacement
%as been set to 0 in the leftact function.
zB1 = gBA .*zA1
zO1 = gOA .*zA1
% b)
zC2 = [-3 ;2 ;pi/8];
%gBC.M(1:2,3) = 0;
zB2 = gBC .* zC2
% ZO2 = gOC * zC2 = (gOB*gBC) * zC2
gOC = gOB*gBC;
%gOC.M(1:2,3) = 0;
zO2 = gOC .* zC2
%% Problem #2
%Given
gBC = SE2([-2;0], -pi/6);
zBB = [2;3;pi/9];
% find zCC
% zCC = Adjoint(gCB,zBB) = gCB*zBB*inv(gCB), but since we are given gBC and
% not gCB we need to write it in terms of gBC which is
% zCC = inv(gBC)*zBB*gBC or Adjoint(inv(gBC),zBB)
zCC = adjoint(inv(gBC),zBB)
%% Problem #3
% Given
gi = [.5,.5,pi/3];
[t,gsol] = ode45(@hw5ode,[0,pi], gi); % for Za = (2,-1,1/3)^T
[t2,gsol2] = ode45(@hw5ode2,[0,pi], gi); % for ZB = (2,-1,0)^T
tt = 1:1:pi;
figure(1)
plot(t,gsol)
hold on
gsol_1 = SE2([gsol(end,1);gsol(end,2)],gsol(end,3));
gsol_2 = SE2([gsol(end-1,1);gsol(end-1,2)],gsol(end-1,3));
plot(gsol_1)
plot(gsol_2)
figure(2)
plot(t2,gsol2)
hold on
gsol2_1 = SE2([gsol2(end,1);gsol2(end,2)],gsol2(end,3));
gsol2_2 = SE2([gsol2(end-1,1);gsol2(end-1,2)],gsol2(end-1,3));
plot(gsol2_1)
plot(gsol2_2)
%% Problem #4
% See Attached SE2 Code
%a)
g = SE2([4;3],3*pi/8);
z = [.7;1.2;.45];
zprime = g.*z
%b)
h = SE2([0,3],-pi/6);
zA = [-.67;3.1;pi/9];
ZB = adjoint(h,zA)
