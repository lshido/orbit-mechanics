% Calculate C for the given position and velocity.
mu_Earth_Moon = 1.2151e-02;

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 1.000];

x = r_vector(1);
y = r_vector(2);
d = sqrt((x+mu_Earth_Moon)^2 + y^2);
r = sqrt((x-1+mu_Earth_Moon)^2 + y^2);
pseudo_U = (1-mu_Earth_Moon)/d + mu_Earth_Moon/r + (x^2+y^2)/2;
v_squared = norm(v_vector)^2;
C = 2*pseudo_U - v_squared;
fprintf("Jacobi Constant C: %d\n", C)