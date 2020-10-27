R_mars = 3397; % km
Gm_mars = 42828.314258067; % km^3/sec^2

TA_maneuver = 150; %[deg]
a_old = 6*R_mars;
e_old = 0.5;
big_omega = 45; % deg
inc = 30; % deg
little_omega = 30; % deg

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

theta = TA_maneuver + little_omega;
fprintf("theta: %.4e deg\n", theta);
fprintf("TA check: %.4e deg\n", theta - little_omega);
delta_v_inertial = [ 0.1 0.25 0.35 ]; % [x^ y^ z^]
fprintf("delta_v_inertial: %.4e km/s x^ %.4e km/s y^ %.4e km/s z^\n", delta_v_inertial);

TA_inertial_to_rotational = [ (cosd(big_omega)*cosd(theta))-(sind(big_omega)*cosd(inc)*sind(theta)) (-cosd(big_omega)*sind(theta))-(sind(big_omega)*cosd(inc)*cosd(theta)) sind(big_omega)*sind(inc);...
                              (sind(big_omega)*cosd(theta))+(cosd(big_omega)*cosd(inc)*sind(theta)) (-sind(big_omega)*sind(theta))+(cosd(big_omega)*cosd(inc)*cosd(theta)) -cosd(big_omega)*sind(inc);...
                              sind(inc)*sind(theta) sind(inc)*cosd(theta) cosd(inc) ];
delta_v_rotational = delta_v_inertial*TA_inertial_to_rotational;
fprintf("delta_v_rotational: %.4e km/s r^ %.4e km/s theta^ %.4e km/s h^\n", delta_v_rotational);

% convert to VNB

TA_rotational_to_vnb = [ sind(FPA_old) cosd(FPA_old) 0 ;...
                         cosd(FPA_old) -sind(FPA_old) 0;...
                         0 0 1 ];
delta_v_vnb = delta_v_rotational*TA_rotational_to_vnb;
fprintf("delta_v_vnb: %.4e km/s v^ %.4e km/s b^ %.4e km/s n^\n", delta_v_vnb(1), delta_v_vnb(2), delta_v_vnb(3));
fprintf("delta_v_vnb magnitude: %.4e km/s\n", norm(delta_v_vnb));                         
fprintf("delta_v_rotational magnitude: %.4e km/s\n", norm(delta_v_rotational));                         
fprintf("delta_v_inertial magnitude: %.4e km/s\n", norm(delta_v_inertial)); 
fprintf("percent of out-of-plane component: %.4e percent\n", delta_v_vnb(3)/norm(delta_v_vnb)*100);
delta_v_r_theta = sqrt((delta_v_rotational(1)^2)+(delta_v_rotational(2)^2));
fprintf("magnitude of projection of delta v on r, theta: %.4e km/s\n", delta_v_r_theta);
beta_angle = atand(delta_v_rotational(3)/delta_v_r_theta);
fprintf("beta_angle: %.4e deg\n", beta_angle);
beta_angle_check = asind(delta_v_rotational(3)/norm(delta_v_rotational));
fprintf("beta_angle_check: %.4e deg\n", beta_angle_check);
beta_angle_check_2 = asind(delta_v_rotational(3)/norm(delta_v_rotational));
fprintf("delta_v_h^: %.4e km/s\n", delta_v_rotational(3));
fprintf("beta_angle_check_2: %.4e deg\n", beta_angle_check_2);
phi_angle = asind(delta_v_rotational(1)/(norm(delta_v_rotational)*cosd(beta_angle)));
fprintf("phi_angle: %.4e deg\n", phi_angle);
phi_angle_check = acosd(delta_v_rotational(2)/(norm(delta_v_rotational)*cosd(beta_angle)));
fprintf("phi_angle_check: %.4e deg\n", phi_angle_check);

% calculate delta_v_v_b
delta_v_v_b = sqrt((delta_v_vnb(1)^2)+(delta_v_vnb(2)^2));
fprintf("magnitude of projection of delta v on v, b: %.4e km/s\n", delta_v_v_b);
alpha_angle = acosd(delta_v_vnb(1)/(norm(delta_v_vnb)*cosd(beta_angle)));
fprintf("alpha_angle: %.4e deg\n", alpha_angle);
alpha_angle_check = asind(delta_v_vnb(2)/(norm(delta_v_vnb)*cosd(beta_angle)));
fprintf("alpha_angle_check: %.4e deg\n", alpha_angle_check);

% express r and v_old in xyz
r_vector_rotational = [ r 0 0 ];
r_vector_inertial = r_vector_rotational*transpose(TA_inertial_to_rotational);
fprintf("r_vector_inertial: %.4e km x^ %.4e km y^ %.4e km z^\n", r_vector_inertial);
fprintf("r magnitude check: %.4e km\n", norm(r_vector_inertial));

v_vector_old_rotational = [ r_dot_old r_theta_dot_old 0 ];
v_vector_old_inertial = v_vector_old_rotational*transpose(TA_inertial_to_rotational);
fprintf("v_vector_old_inertial: %.4e km/s x^ %.4e km/s y^ %.4e km/s z^\n", v_vector_old_inertial);
fprintf("v old magnitude check: %.4e km\n", norm(v_vector_old_inertial));

v_vector_new_inertial = v_vector_old_inertial + delta_v_inertial;
fprintf("v_vector_new_inertial: %.4e km/s x^ %.4e km/s y^ %.4e km/s z^\n", v_vector_new_inertial);






