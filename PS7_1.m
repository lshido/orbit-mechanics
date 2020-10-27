R_mars = 3397; % km
Gm_mars = 42828.314258067; % km^3/sec^2

TA_maneuver = 150; %[deg]
a_old = 6*R_mars;
e_old = 0.5;

fprintf("Semi major axis before maneuver: %.4e km\n", a_old);

specific_energy_old = -Gm_mars/(2*a_old);
fprintf("specific energy before maneuver: %.4e km^2/s^2\n", specific_energy_old);

p_old = a_old*(1-e_old^2);
fprintf("semilatus rectum before maneuver: %.4e km\n", p_old);

r = p_old/(1+e_old*cosd(TA_maneuver));
fprintf("radius at maneuver: %.4e km\n", r);

% calc v
v_old = sqrt(2*(specific_energy_old+(Gm_mars/r)));
fprintf("Velocity mag before maneuver: %.4e km/sec\n", v_old);

% calc r_dot 
r_dot_old = sqrt((v_old^2)-(p_old*Gm_mars/(r^2)));
fprintf("r_dot mag before maneuver: %.4e km/sec\n", r_dot_old);

% calc r_theta_dot
r_theta_dot_old = sqrt(p_old*Gm_mars)/r;
fprintf("r_theta_dot mag before maneuver: %.4e km/sec\n", r_theta_dot_old);

% calc FPA
FPA_old = atand(r_dot_old/r_theta_dot_old); % positive because ascending
fprintf("FPA before maneuver: %.4e deg\n", FPA_old);

% convert to rotational
big_omega = 45; % deg
inc = 30; % deg
little_omega = 30; % deg
theta = TA_maneuver + little_omega;
delta_v_inertial = [ 0.1 0.25 0.35 ]; % [x^ y^ z^]

TA_inertial_to_rotational = [ (cosd(big_omega)*cosd(theta))-(sind(big_omega)*cosd(inc)*sind(theta)) (-cosd(big_omega)*sind(theta))-(sind(big_omega)*cosd(inc)*cosd(theta)) sind(big_omega)*sind(i);...
                              (sind(big_omega)*cosd(theta))+(cosd(big_omega)*cosd(inc)*sind(theta)) (-sind(big_omega)*sind(theta))+(cosd(big_omega)*cosd(inc)*cosd(theta)) -cosd(big_omega)*sind(inc);...
                              sind(inc)*sind(theta) sind(inc)*cosd(theta) cosd(i) ];
delta_v_rotational = delta_v_inertial*TA_inertial_to_rotational;
fprintf("delta_v_rotational: %.4e km/s r^ %.4e km/s theta^ %.4e km/s h^\n", delta_v_rotational(1), delta_v_rotational(2), delta_v_rotational(3));

% convert to VNB

TA_rotational_to_vnb = [ sind(FPA_old) cosd(FPA_old) 0 ;...
                         cosd(FPA_old) -sind(FPA_old) 0;...
                         0 0 1 ];
delta_v_vnb = delta_v_rotational*TA_rotational_to_vnb;
fprintf("delta_v_vnb: %.4e km/s v^ %.4e km/s b^ %.4e km/s n^\n", delta_v_vnb(1), delta_v_vnb(2), delta_v_vnb(3));
                         

