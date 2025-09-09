R_earth = 6378.1363; % [km]

% Solve for a and e
r_p = 1.5*R_earth;
r_a = 200*R_earth;
fprintf("radius at periapsis: %.4e km\n", r_p);
fprintf("radius at apoapsis: %.4e km\n", r_a);

syms e a
equations = [a*(1-e) == r_p, a*(1+e) == r_a];
variables = [a e];
[solved_a, solved_e] = solve(equations, variables);
fprintf("semi-major axis: %.4e km\n", solved_a);
fprintf("eccentricity: %.4e\n", solved_e);

% ------------------------------------------------

r_sun_earth = 149597898;
r_earth_moon = -384400;
r_earth_sc = -r_a;

% known Gm values
Gm_earth = 398600.4415;
Gm_sun = 132712440017.99;
Gm_moon = 4902.8005821478;

r_earth_sun = -r_sun_earth;
r_sc_sun = -r_sun_earth - r_earth_sc;
r_sc_moon = r_earth_moon - r_earth_sc;

dominant = -(Gm_earth*r_earth_sc)/((abs(r_earth_sc))^3);
fprintf('dominant = %.5e km/s^2\n\n', dominant);

% due to earth
direct_sun = Gm_sun*r_sc_sun/(abs(r_sc_sun))^3;
indirect_sun = Gm_sun*r_earth_sun/(abs(r_earth_sun))^3;
net_pert_sun = direct_sun - indirect_sun;
fprintf('direct_sun = %.5e km/s^2\n', direct_sun);
fprintf('indirect_sun = %.5e km/s^2\n', indirect_sun);
fprintf('net_pert_sun = %.5e km/s^2\n\n', net_pert_sun);

% due to sun
direct_moon = Gm_moon*r_sc_moon/(abs(r_sc_moon))^3;
indirect_moon = Gm_moon*r_earth_moon/(abs(r_earth_moon))^3;
net_pert_moon = direct_moon - indirect_moon;
fprintf('direct_moon = %.5e km/s^2\n', direct_moon);
fprintf('indirect_moon = %.5e km/s^2\n', indirect_moon);
fprintf('net_pert_moon = %.5e km/s^2\n\n', net_pert_moon);