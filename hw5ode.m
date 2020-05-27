function gdot = hw5ode(t, g)
xi = [2; -1; 1/3];
th = g(3); 
gdot = [ cos(th) , -sin(th) , 0 ; sin(th), cos(th), 0; 0, 0, 1 ]*xi;
end 