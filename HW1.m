%% ECE 4560 HW#1
% Submitted By: Jesse Austin Stringfellow
% Due Wednesday, Sept 4, 2019
%% Problem 1
p0 = [6;0]; %global var for problem 1
%a)
T = [10;25]; %Translation
theta = (pi/6); %Rotation
rotMat = rotMatrix(theta); %Function Makes Rotation Matrix from theta
p_a = rotMat*p0 + T %Rotates and then translates a point
%b)
Tb = [21;7];
thetab = pi;
rotMatb = rotMatrix(thetab);
p_b = rotMatb*p0 + Tb
%% Problem 2
% See at back of Homework
%% Problem 3
q = [1;0];

gAO = [5,12,pi/3];
d = [gAO(1);gAO(2)];
R = rotMatrix(gAO(3));
gAO_in_dR = [d,R]

gBO = [2,-1,pi];
d = [gBO(1);gBO(2)];
R = rotMatrix(gBO(3));
gBO_in_dR = [d,R]

%% Problem 4

qAO = gAO_in_dR(1:2,2:3)*q + gAO_in_dR(1:2,1)

%% Problem 5
qBO = gBO_in_dR(1:2,2:3)*q + gBO_in_dR(1:2,1)
%% Problem 6
% gCO is found by multiplying gAO * gCA as per class
gCA = [-2,-2,pi/6];
d = [gCA(1);gCA(2)];
R = rotMatrix(gCA(3));
gCA_in_dR = [d,R]

gCO_in_dR = [gAO_in_dR(1:2,1) + gAO_in_dR(1:2,2:3)*gCA_in_dR(1:2,1),gAO_in_dR(1:2,2:3)*gCA_in_dR(1:2,2:3)]
qCO = gCO_in_dR(1:2,2:3)*q + gCO_in_dR(1:2,1)

%% Problem 7
qnew = [.6428; -.7660];
qAO = gAO_in_dR(1:2,2:3)*qnew + gAO_in_dR(1:2,1)
