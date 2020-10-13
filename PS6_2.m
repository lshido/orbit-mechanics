R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

a = 4*R_earth;
fprintf("Semi major axis: %.4e km\n", a);
e = 0.4;
fprintf("eccentricity: %.4e\n", e);
TA = 135; %[deg]


% calc p
p = a*(1-e^2);
fprintf("semilatus rectum: %.4e km\n", p);

% calc r at maneuver point
r = p/(1+e*cosd(TA));
fprintf("r at maneuver point: %.4e km\n", r);

% calc v new
v_old = sqrt(2*((Gm_earth/r) - (Gm_earth/(2*a))));
fprintf("velocity magnitude before maneuver: %.4e km/sec\n", v_old);

% calc gamma_old
gamma_old = acosd(sqrt(p*Gm_earth)/(r*v_old));
fprintf("FPA before maneuver: %.4e deg\n", gamma_old);

% calc v vector old
v_vector_old = [v_old*sind(gamma_old) v_old*cosd(gamma_old)];
fprintf("Velocity vector before maneuver: %.4e km/s %.4e km/s\n", v_vector_old(1), v_vector_old(2));