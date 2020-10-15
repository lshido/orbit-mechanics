R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

a = 4*R_earth;
fprintf("Semi major axis: %.4e km\n", a);
e = 0.4;
fprintf("eccentricity: %.4e\n", e);
TA = 135; %[deg]

d_v = 0.90; % [km/sec]
alpha_angle = 45; % [deg]

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

% express r and v in inertial frame
TM_inertial = [ cosd(TA) sind(TA); -sind(TA) cosd(TA)];
r_inertial = [r 0] * TM_inertial;
v_inertial_old = v_vector_old * TM_inertial;
fprintf("r_inertial: %.4e km e^ %.4e km p^\n", r_inertial(1), r_inertial(2));
fprintf("v_inertial before maneuver: %.4e km/s e^ %.4e km/s p^\n", v_inertial_old(1), v_inertial_old(2));

% express r and v in VNB frame
TM_vnb = [ cosd(90-gamma_old) sind(90-gamma_old); cosd(gamma_old) -sind(gamma_old)];
r_vnb = [r 0] * TM_vnb;
v_vnb_old = v_vector_old * TM_vnb;
fprintf("r_vnb: %.4e km v^ %.4e km b^\n", r_vnb(1), r_vnb(2));
fprintf("v_vnb before maneuver: %.4e km/s v^ %.4e km/s b^\n", v_vnb_old(1), v_vnb_old(2));

% find v_new
v_new = sqrt((v_old^2) + (d_v^2) - (2*v_old*d_v*cosd(180-alpha_angle)));
fprintf("v_new magnitude: %.4e km/sec\n", v_new);

% find d_gamma
d_gamma = acosd((d_v^2-v_old^2-v_new^2)/(-2*v_old*v_new));
fprintf("d_gamma: %.4e deg\n", d_gamma);
gamma_new = gamma_old + d_gamma;
fprintf("gamma_new: %.4e deg\n", gamma_new);

% find vector new
v_vector_new = [v_new*sind(gamma_new) v_new*cosd(gamma_new)];
fprintf("Velocity vector after maneuver: %.4e km/s r^ %.4e km/s theta^\n", v_vector_new(1), v_vector_new(2));

% express r and v in inertial frame
v_inertial_new = v_vector_new * TM_inertial;
fprintf("v_inertial after maneuver: %.4e km/s e^ %.4e km/s p^\n", v_inertial_new(1), v_inertial_new(2));
v_vnb_new = v_vector_new * TM_vnb;
fprintf("v_vnb after maneuver: %.4e km/s v^ %.4e km/s b^\n", v_vnb_new(1), v_vnb_new(2));

% calc h (specific angular momentum)
h = r*v_vector_new(2);
fprintf("specific angular momentum: %.4e km^2/sec\n", h);
p_new = h^2/Gm_earth;
fprintf("semilatus rectum after maneuver: %.4e km\n", p_new);

% calc new specific energy
specific_energy_new = (v_new^2/2) - (Gm_earth/r);
fprintf("specific energy after maneuver: %.4e km^2/sec^2\n", specific_energy_new);

% calc new semi major axis
a_new = -Gm_earth/(2*specific_energy_new);
fprintf("semi major axis after maneuver: %.4e km\n", a_new);

% calc new eccentricity
e_new = sqrt(1-(p_new/a_new));
fprintf("eccentricity after maneuver: %.4e\n", e_new);

% calc new period
period_new = 2*pi*sqrt((a_new^3)/Gm_earth);
fprintf("period after maneuver: %.4e sec\n", period_new);
fprintf("period after maneuver: %.4e hours\n", period_new/3600);
fprintf("period after maneuver: %.4e days\n", period_new/3600/24);

% calc new periapsis radius
rp_new = a_new*(1-e_new);
fprintf("rp after maneuver: %.4e km\n", rp_new);

% calc EA
EA_new = 2*atan(sqrt((1-e_new)/(1+e_new))*tand(TA/2));
fprintf("Eccentric Anomaly: %.4e rad\n", EA_new);
fprintf("Eccentric Anomaly: %.4e deg\n", rad2deg(EA_new));

% calc MA
MA_new = EA_new-e_new*sin(EA_new);
fprintf("Mean Anomaly: %.4e rad\n", MA_new);
fprintf("Mean Anomaly: %.4e deg\n", rad2deg(MA_new));

%calc time since periapsis
t_since_tp = MA_new*sqrt((a_new^3)/Gm_earth);
fprintf("time since periapsis: %.4e sec\n", t_since_tp);
fprintf("time since periapsis: %.4e hours\n", t_since_tp/3600);
fprintf("time since periapsis: %.4e days\n", t_since_tp/3600/24);



