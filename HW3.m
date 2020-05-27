%% ECE 4560 HW #3 
% Submitted by Jesse Austin Stringfellow 

% Due Sept. 18, 2019
%% Problem 1
% Transforms from HW 1
gOO = SE2([0;0],0);
gOA = SE2([5;12], pi/3); 
gOB = SE2([2;-1], pi);
gAC = SE2([-2;-2], pi/6);
gOC = gOA*gAC;
% Points from HW 1 
qOA = [5.5; 12.866]; 
qOB = [1;-1]; 
qOC = [5.7321;10.2679]; 
qOAm = [5.9848;12.1737];

figure(1); clf; 
gOO.plot('O','r'); 
hold on; %Used to allow multiple plots
gOA.plot('A');
gOB.plot('B'); 
gOC.plot('C','m');

%plots points with the nudged hammer location in Cyan
plot(qOA(1), qOA(2), 'b*', 'MarkerSize', 8);
plot(qOB(1), qOB(2), 'b*', 'MarkerSize', 8);
plot(qOC(1), qOC(2), 'b*', 'MarkerSize', 8); 
plot(qOAm(1), qOAm(2), 'c*', 'MarkerSize', 8);
hold off;
axis equal;
axis([-5 15 -5 15]);

%% Problem 2
% View back of HW for SE2 Code
%g1 and g2 are given
g1 = SE2([1;2], pi/3);
g2 = SE2([-2;1], pi/6);
%verification Code
inverse_g1 = inv(g1) 
g3 = g1*g2
%% Problem 3
% First two parts are attached to the back in handwriting.
%Part a
    %See attached notepaper
%Part b
    %See attached notepaper
%Part c
%Other way is to plug values into the vector form in part a)
gm1 = SE2([0;0], pi); 
gm2 = SE2([1.5;0], pi/8); 
gm3 = SE2([0.5;0], -pi/4);
gm4 = SE2([0.3;0], pi/8);
ge = gm1*gm2*gm3*gm4
%Part d
gma= gm1*gm2; % Corresponds to where the second joint is in relation to the origin
gmb = gm1*gm2*gm3; % Corresponds to where the third joint is in relation to the origin
figure(2); clf; 
gm1.plot('1','r'); % RED Plots first joint
hold on; 
gma.plot('2','b'); % BLUE Plots second joint with relation to origin
gmb.plot('3','g'); % GREEN Plots third joint with relation to origin
ge.plot('4','c'); %  CYAN Plots end-effector location with relation to origin
hold off;
axis equal;
axis([-4 3 -4 3]);



