%% ECE 4560 HW #2
% Submitted by: Jesse Austin Stringfellow

%% Problem 7
%% a)
tspan = [0, 20];
a0 = [-3*pi/2; pi/6];
[tsol, a] = ode45( @hw2, tspan, a0);

figure(1);
plot(tsol, a(:,1));
hold on
plot(tsol,a(:,2));
xlabel('t');
ylabel('x');
grid on;
hold off


%% b)
%Animation works, not sure how to show that when publishing. Code was used
%and adapted from Dr. Vela's webpage.
l1 = [1;1/2];
ti = 0;
tf = 20;
nframes = 100;
tvect = linspace(ti,tf,nframes);
figure(2);
for tT = tvect
  alphaT = transpose(interp1(tsol, a, tT));
  planarR2_display(alphaT, l1);
  drawnow;
end

%% Functions
%Function used in @hw2 for Part a
function adot = f(t, a)

adot = zeros(2,1);
adot(1) = (1/3)*cos(t); 
adot(2) = -(1/4)*sin(t);
end
