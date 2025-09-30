% spacecraft in the Earth-Moon system
% CR3BP

% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% position and velocity in NON-DIMENSIONAL units
r_vector = [-0.270 -0.420];
v_vector = [0.300 -1.000];

r_dim_vector = r_vector * l_char *1000;
fprintf("dimensional position: %7e m\n", r_dim_vector)

v_dim_vector = v_vector * l_char/t_char*1000;
fprintf("dimensional vector: %7e m/s\n", v_dim_vector)


