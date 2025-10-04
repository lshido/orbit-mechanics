clear all

mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("mu %d\n", mu)

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_x = 0.3;
x = r_vector(1);
y = r_vector(2);
C = 3.150;

d = sqrt((x+mu).^2 + y.^2);
r = sqrt((x-1+mu).^2 + y.^2);
x_y_sq = (x.^2+y.^2)./2;
term_1 = (1-mu)./d;
term_2 = mu./r;
pseudo_U = term_1 + term_2 + x_y_sq;
% v_squared = w(:,3).^2 + w(:,4).^2;
v_squared = 2*pseudo_U - C;
v_y = sqrt(v_squared - v_x^2);


v_squared_check = v_x^2 + v_y^2;
C_check = 2*pseudo_U - v_squared;

fprintf("velocity_x: %d, velocity_y: %d\n", v_x, v_y)
fprintf("JC Check: %d\n", C_check)
