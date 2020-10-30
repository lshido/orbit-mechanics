R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

%----------Single Maneuver-----------

% Current orbit
r = 200 + R_earth;
e = 0;
a = r;
p = r;

% Conditions prior to maneuver
fprintf("r: %.4e km\n", r);
specific_energy_old = -Gm_earth/(2*a);
fprintf("specific_energy_old: %.4e km^2/s^2\n", specific_energy_old);
v_old = sqrt(2*(specific_energy_old+(Gm_earth/r)));
fprintf("v_old: %.4e km/s\n", v_old);

% Conditions after maneuver
inc = 57;
v_new = v_old;
delta_v = 2*v_new*sind(inc/2);
fprintf("delta_v: %.4e km/s\n", delta_v);

% calculate beta
psi_angle = (180-inc)/2;
fprintf("psi_angle: %.4e deg\n", psi_angle);
beta_angle = 180-psi_angle;
fprintf("beta_angle: %.4e deg\n", beta_angle);

%--------------Bi-elliptic transfer--------------
r_2 = 55*R_earth;

% determine transfer ellipse
a_transfer_1 = 1/2*(r + r_2);
fprintf("a_transfer_1: %.4e km\n", a_transfer_1);
e_transfer_1 = 1 - (r/a_transfer_1);
fprintf("e_transfer_1: %.4e \n", e_transfer_1);

% conditions after maneuver at 1
v_1_after = sqrt(2*((Gm_earth/r)-(Gm_earth/(2*a_transfer_1))));
fprintf("v_1_after: %.4e km/sec\n", v_1_after);

% calculate delta v 1
delta_v_1 = v_1_after - v_old;
fprintf("delta_v_1: %.4e km/sec\n", delta_v_1);

% conditions before maneuver at 2
v_2_before = sqrt(2*((Gm_earth/r_2)-(Gm_earth/(2*a_transfer_1))));
fprintf("v_2_before: %.4e km/sec\n", v_2_before);

% Conditions after maneuver
inc = 57;
delta_v_2 = 2*v_2_before*sind(inc/2);
fprintf("delta_v_2: %.4e km/s\n", delta_v_2);

% calculate beta
psi_angle_2 = (180-inc)/2;
fprintf("psi_angle_2: %.4e deg\n", psi_angle_2);
beta_angle_2 = 180-psi_angle;
fprintf("beta_angle_2: %.4e deg\n", beta_angle_2);
delta_v_2_vector = [ delta_v_2*cosd(beta_angle) delta_v_2*sind(beta_angle) 0 ];
fprintf("delta_v_2_vector: %.4e km/s V^ %.4e km/s N^ %.4e km/s B^\n", delta_v_2_vector);


% conditions before maneuver at 3
a_transfer_2 = a_transfer_1;
v_3_before = sqrt(2*((Gm_earth/r)-(Gm_earth/(2*a_transfer_2))));
fprintf("v_3_before: %.4e km/sec\n", v_3_before);

% conditions required after maneuver in final orbit
v_3_after = sqrt(Gm_earth/r);
fprintf("v_3_after: %.4e km/sec\n", v_3_after);
delta_v_3 = v_3_after - v_3_before;
fprintf("delta_v_3: %.4e km/sec\n", delta_v_3);

% calculate total delta_v
delta_v_total = abs(delta_v_1) + abs(delta_v_2) + abs(delta_v_3);
fprintf("delta_v_total: %.4e km/sec\n", delta_v_total);

% calculate TOF
TOF_bielliptic = 2*pi*sqrt(a_transfer_1^3/Gm_earth);
fprintf("TOF_bielliptic in sec: %.8e sec\n", TOF_bielliptic);
fprintf("TOF_bielliptic in days: %.8e days\n", TOF_bielliptic/3600/24);
fprintf("TOF_bielliptic in years: %.8e years\n", TOF_bielliptic/3600/24/365);

% delta v cost savings
delta_v_savings = delta_v - delta_v_total;
fprintf("delta_v_savings: %.4e km/sec\n", delta_v_savings);
