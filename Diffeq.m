
%define our differential funtion
function alphadot = Diffeq(t, alpha)

vdes = [-0.5;0.5];
J = pJac(alpha);
 alphadot = pinv(J)*vdes;
end
