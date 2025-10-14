clear all

%================================DEFINE FUNCTIONS=====================
%============================END DEFINE FUNCTIONS=====================

%====================DEFINE SYSTEM====================================
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("mu %d\n", mu)

% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;

% Calculate characteristic time
t_char = sqrt(l_char^3/(mu_Earth+mu_Moon));

% position and velocity in NON-DIMENSIONAL units
r_vector = [0.488 0.200];
v_vector = [-0.880 0.200];

%====================DEFINE SYSTEM====================================

%====================PRINT IMPORTANT NUMBERS==========================
fprintf("characteristic time: %f sec\n", t_char)
fprintf("characteristic length: %f sec\n", l_char)
fprintf("Initial x: %f km\n", r_vector(1)*l_char)
fprintf("Initial y: %f km\n", r_vector(2)*l_char)
fprintf("Initial v_x: %f m/s\n", v_vector(1)*l_char/t_char*1000)
fprintf("Initial v_y: %f m/s\n", v_vector(2)*l_char/t_char*1000)

%================END PRINT IMPORTANT NUMBERS==========================