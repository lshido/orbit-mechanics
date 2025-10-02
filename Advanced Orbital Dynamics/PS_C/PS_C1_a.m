% Calculate C for the given position and velocity.
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);

% mu_Earth_Moon = 1.2151e-02;

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 -1.000];

x = r_vector(1);
y = r_vector(2);
d = sqrt((x+mu)^2 + y^2);
r = sqrt((x-1+mu)^2 + y^2);
pseudo_U = (1-mu)/d + mu/r + (x^2+y^2)/2;
v_squared = norm(v_vector)^2;
C = 2*pseudo_U - v_squared;
fprintf("Jacobi Constant C: %d\n", C)

fprintf("Mu: %d\n", mu)

term1_num = 2*(1-mu);
fprintf("First Term Numerator: %d\n", term1_num)

term2_num = 2*mu;
fprintf("Second Term Numerator: %d\n", term2_num)

l = 1+mu;
fprintf("1+mu: %d\n", term2_num)


 