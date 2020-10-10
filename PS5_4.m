R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

r_inertial_1 = [0.15*R_earth -1.44*R_earth -0.65*R_earth];
v_inertial_1 = [6.62 2.7 -1.56];

r_1 = norm(r_inertial_1);
v_1 = norm(v_inertial_1);
fprintf("r at t1: %.4e km\n", r_1);
fprintf("v at t1: %.4e km/sec\n", v_1);

% specific energy
energy = ((v_1^2)/2) - (Gm_earth/r_1);
fprintf("Specific Energy: %.4e km^2/sec^2\n", energy);

% Semi major axis
a = - Gm_earth/(2*energy);
fprintf("semi major axis: %.4e km\n", a);

% escape speed
escape = sqrt((2*Gm_earth)/r_1);
fprintf("escape speed: %.4e km/sec\n", escape);

% h
h_vector = cross(r_inertial_1, v_inertial_1);
h = norm(h_vector);
fprintf("specific angular momentum: %.4e km^2/sec\n", h);

% h unit vector
h_unit = h_vector/h;
fprintf("h unit vector: %.4e x^ %.4e y^ %.4e z^ km^2/sec\n", h_unit(1), h_unit(2), h_unit(3));
fprintf("h unit vector check: %.4e km^2/sec\n", norm(h_unit));

% semilatus rectum
p = h^2/Gm_earth;
fprintf("semlatus rectum: %.4e km\n", p);

% eccentricity
e = sqrt(1-(p/a));
fprintf("eccentricity: %.4e \n", e);

% true anomaly
TA_1 = acosd((1/e)*((p/r_1)-1));
fprintf("TA at t1: %.4e degrees\n", TA_1);

% inclination
inclination = acosd(h_unit(3));
fprintf("inclination: %.4e degrees\n", inclination);

% big omega
big_omega_x_value = h_unit(1)/sind(inclination);
big_omega_x = asind(big_omega_x_value);
fprintf("big omega x: %.4e degrees\n", big_omega_x);
fprintf("big omega x value: %.4e\n", big_omega_x_value);

% big omega
big_omega_y_value = -h_unit(2)/sind(inclination);
big_omega_y = acosd(big_omega_y_value);
fprintf("big omega y: %.4e degrees\n", big_omega_y);
fprintf("big omega y value: %.4e\n", big_omega_y_value);

% r unit vector
r_unit = r_inertial_1/r_1;
fprintf("r unit vector: %.4e x^ %.4e y^ %.4e z^ km\n", r_unit(1), r_unit(2), r_unit(3));
fprintf("r unit vector check: %.4e km\n", norm(r_unit));

% theta unit vector
theta_unit = cross(h_unit, r_unit);
fprintf("theta unit vector: %.4e x^ %.4e y^ %.4e z^ km\n", theta_unit(1), theta_unit(2), theta_unit(3));
fprintf("theta unit vector check: %.4e km\n", norm(theta_unit));

% theta
theta_r_value = r_unit(3)/sind(inclination);
theta_r = asind(theta_r_value);
fprintf("theta r: %.4e degrees\n", theta_r);
fprintf("theta r value: %.4e\n", theta_r_value);

% theta
theta_t_value = theta_unit(3)/sind(inclination);
theta_t = acosd(theta_t_value);
fprintf("theta t: %.4e degrees\n", theta_t);
fprintf("theta t value: %.4e\n", theta_t_value);

% velocity r component
r_dot_1 = dot(v_inertial_1, r_unit);
fprintf("r_dot_1: %.4e km/sec\n", r_dot_1);

r_theta_dot_1 = dot(v_inertial_1, theta_unit);
fprintf("r_theta_dot_1: %.4e km/sec\n", r_theta_dot_1);

% FPA1
FPA_1 = atand(r_dot_1/r_theta_dot_1);
fprintf("FPA_1: %.4e degrees\n", FPA_1);

% calc EA1
EA_1 = 2*atan(sqrt((1-e)/(1+e))*tand(-TA_1/2));
fprintf("Eccentric Anomaly at t1: %.4e rad\n", EA_1);
fprintf("Eccentric Anomaly at t1: %.4e deg\n", rad2deg(EA_1));

% calc MA1
MA_1 = EA_1 - e*sin(EA_1);
fprintf("Mean Anomaly at t1: %.4e rad\n", MA_1);
fprintf("Mean Anomaly at t1: %.4e deg\n", rad2deg(MA_1));

%calc time since periapsis
t1_since_tp = MA_1*sqrt((a^3)/Gm_earth);
fprintf("time since periapsis: %.4e sec\n", t1_since_tp);
fprintf("time since periapsis: %.4e hours\n", t1_since_tp/3600);
fprintf("time since periapsis: %.4e days\n", t1_since_tp/3600/24);

%calc period
period = 2*pi*sqrt((a^3)/Gm_earth);
fprintf("period: %.4e sec\n", period);
fprintf("period: %.4e hours\n", period/3600);
fprintf("period: %.4e days\n", period/3600/24);







