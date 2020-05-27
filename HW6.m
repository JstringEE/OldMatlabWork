%% HW #6
% Submitted by Jesse Austin Stringfellow,    
% Due Oct. 30, 2019
%% Problem #1
g1 = SE3 ([1;2;3] ,[1 0 0;0 0 -1; 0 1 0]);
g2 = SE3 ([2; -1; -1] , [sqrt(2)/2 0 -sqrt(2)/2;0 1 0; sqrt(2)/2 0 sqrt(2)/2]);
p1 = [4;5;6];
p2 = [4;5;6;1];
g3 = g1 * g2
g4 = g3 .* p1 
g5 = g3 .* p2
invg3 = inv(g3)
%% Problem #2 See attached page
