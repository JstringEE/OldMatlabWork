%% HW #7
% Submitted by Jesse Austin Stringfellow,    
% Due Nov. 13, 2019
%% Problem #1
planar_r3(1.5,1,.3)
%% Problem #2
%% a)
syms a1 a2 a3 a4 l0 l1 l2 l3
g1 = SE3([0;0;l0],[cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1])
g2 = SE3([0;0;0],[1 0 0; 0 cos(a2) -sin(a2); 0 sin(a2) cos(a2)])
g3 = SE3([0;l1;0],[1 0 0; 0 cos(a3) -sin(a3); 0 sin(a3) cos(a3)])
g4 = SE3([0;l2;0],[1 0 0; 0 cos(a4) -sin(a4); 0 sin(a4) cos(a4)])
g5 = SE3([0;l3;0],eye(3))
ge = g1*g2*g3*g4*g5
getRotationMatrix(ge)
getTranslation(ge)
%% b)
reachpartc(1,3/4,1/2) %Modified planar_r3 that reflects the
                      %joint limits and the translation matrix
%% c) 
a1 = pi/3;a2=pi/3;a3=-pi/4;
l0 = 1; l1 = 3/4; l2 = 1/2;
g1 = SE3([0;0;l0],[cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1]);
g2 = SE3([0;0;0],[1 0 0; 0 cos(a2) -sin(a2); 0 sin(a2) cos(a2)]);
g3 = SE3([0;l1;0],[cos(a3) -sin(a3) 0; sin(a3) cos(a3) 0; 0 0 1]);
g4 = SE3([0;l2;0],eye(3));
ge = g1*g2*g3*g4