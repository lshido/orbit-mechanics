R_jupiter = 71492;
Gm_jupiter = 126712767.8578; %km^3/sec^2


fprintf("-----------Problem 3a: Remus Orbit---------\n")
a_remus = 15*R_jupiter;
e_remus = 0.25;
R_remus = 3000;
Gm_remus = 1e5;
rp_remus = a_remus*(1-e_remus);
ra_remus = a_remus*(1+e_remus);
p_remus = a_remus*(1-e_remus^2);
fprintf("a_remus: %.4e km\n", a_remus);
fprintf("e_remus: %.4e\n", e_remus);
fprintf("p_remus: %.4e km\n", p_remus);
fprintf("rp_remus: %.4e km\n", rp_remus);
fprintf("ra_remus: %.4e km\n", ra_remus);

fprintf("-----------Problem 3a: Spacecraft Orbit---------\n")
e_spacecraft = 0.5;
rp_spacecraft = 7.5*R_jupiter;
a_spacecraft = rp_spacecraft/(1-e_spacecraft);
ra_spacecraft = a_spacecraft*(1+e_spacecraft);
p_spacecraft = a_spacecraft*(1-e_spacecraft^2);
fprintf("a_spacecraft: %.4e km\n", a_spacecraft);
fprintf("e_spacecraft: %.4e\n", e_spacecraft);
fprintf("p_spacecraft: %.4e km\n", p_spacecraft);
fprintf("rp_spacecraft: %.4e km\n", rp_spacecraft);
fprintf("ra_spacecraft: %.4e km\n", ra_spacecraft);

fprintf("-----------Problem 3a: Encounter location---------\n")
r_encounter = a_spacecraft; % at end of spacecraft minor axis
TA_spacecraft = acosd(-e_spacecraft);
TA_remus = acosd((p_remus/(r_encounter*e_remus))-(1/e_remus));
delta_small_omega_remus_spacecraft = -TA_remus + TA_spacecraft;
fprintf("r_encounter: %.4e km\n", r_encounter);
fprintf("TA_spacecraft: %.4e deg\n", TA_spacecraft);
fprintf("TA_remus: %.4e deg\n", TA_remus);
fprintf("delta_small_omega_remus_spacecraft: %.4e deg\n", delta_small_omega_remus_spacecraft);

fprintf("-----------Problem 3b: Spacecraft at Encounter---------\n")
energy_encounter_sc = -Gm_jupiter/(2*a_spacecraft);
v_encounter_old_mag_sc = sqrt(2*(energy_encounter_sc + (Gm_jupiter/r_encounter)));
FPA_encounter_old_sc = acosd(sqrt(p_spacecraft*Gm_jupiter)/(r_encounter*v_encounter_old_mag_sc));
v_encounter_old_rhat = sqrt((v_encounter_old_mag_sc^2) - ((p_spacecraft*Gm_jupiter)/r_encounter^2));
v_encounter_old_thetahat = sqrt(p_spacecraft*Gm_jupiter)/r_encounter;
FPA_encounter_old_sc_check = atand(v_encounter_old_rhat/v_encounter_old_thetahat);
v_encounter_old_sc = [ v_encounter_old_mag_sc*sind(FPA_encounter_old_sc) v_encounter_old_mag_sc*cosd(FPA_encounter_old_sc) ];

fprintf("energy_encounter_sc: %.4e km^2/s^2\n", energy_encounter_sc);
fprintf("v_encounter_old_mag_sc: %.4e km/s\n", v_encounter_old_mag_sc);
fprintf("FPA_encounter_old_sc: %.4e deg\n", FPA_encounter_old_sc);
fprintf("FPA_encounter_old_sc_check: %.4e deg\n", FPA_encounter_old_sc_check);
fprintf("v_encounter_old_sc: %.4e km/s r^ %.4e km/s theta^\n", v_encounter_old_sc);

plot_eph(e_remus, p_remus, delta_small_omega_remus_spacecraft, 0, 360)
plot_eph(e_spacecraft, p_spacecraft, 0, 0, 360)

fprintf("-----------Problem 3b: Remus at Encounter---------\n")
energy_encounter_remus = -Gm_jupiter/(2*a_remus);
v_encounter_old_mag_remus = sqrt(2*(energy_encounter_remus + (Gm_jupiter/r_encounter)));
FPA_encounter_old_remus = acosd(sqrt(p_remus*Gm_jupiter)/(r_encounter*v_encounter_old_mag_remus));
v_encounter_old_remus = [ v_encounter_old_mag_remus*sind(FPA_encounter_old_remus) v_encounter_old_mag_remus*cosd(FPA_encounter_old_remus) ];

fprintf("energy_encounter_remus: %.4e km^2/s^2\n", energy_encounter_remus);
fprintf("v_encounter_old_mag_remus: %.4e km/s\n", v_encounter_old_mag_remus);
fprintf("FPA_encounter_old_remus: %.4e deg\n", FPA_encounter_old_remus);
fprintf("v_encounter_old_remus: %.4e km/s r^ %.4e km/s theta^\n", v_encounter_old_remus);

fprintf("-----------Problem 3d: Compute post-encounter---------\n")
rp_hyp = R_remus + 1500;
v_inf_arr = v_encounter_old_sc - v_encounter_old_remus;
v_inf_arr_mag = norm(v_inf_arr);
a_hyp = Gm_remus/(v_inf_arr_mag^2);
e_hyp = (rp_hyp/a_hyp) + 1;
turn_angle_hyp = asind(1/e_hyp);
FBA_hyp = 2*turn_angle_hyp;
fprintf("rp_hyp: %.4e km\n", rp_hyp);
fprintf("v_inf_arr: %.4e km/s r^ %.4e km/s theta^\n", v_inf_arr);
fprintf("v_inf_arr_mag: %.4e km/s\n", v_inf_arr_mag);
fprintf("a_hyp: %.4e km\n", a_hyp);
fprintf("e_hyp: %.4e \n", e_hyp);
fprintf("turn_angle_hyp: %.4e deg\n", turn_angle_hyp);
fprintf("FBA_hyp: %.4e deg\n", FBA_hyp);

fprintf("-----------Problem 3d: Compute post-encounter---------\n")
delta_v_eq_mag = 2*v_inf_arr_mag*sind(FBA_hyp/2);
nu_angle = FPA_encounter_old_sc - FPA_encounter_old_remus;
rho_angle = asind((v_encounter_old_mag_sc*sind(nu_angle))/v_inf_arr_mag);
rho_angle_check = acosd(((v_encounter_old_mag_sc^2) - (v_encounter_old_mag_remus^2) - (v_inf_arr_mag^2))/(-2*v_encounter_old_mag_remus*v_inf_arr_mag));
sigma_angle = FBA_hyp - rho_angle;
v_encounter_new_mag_sc = sqrt((v_encounter_old_mag_remus^2)+(v_inf_arr_mag^2)-(2*v_encounter_old_mag_remus*v_inf_arr_mag*cosd(sigma_angle)));
kappa_angle = asind((v_inf_arr_mag*sind(sigma_angle))/v_encounter_new_mag_sc);
kappa_angle_check = acosd(((v_inf_arr_mag^2) - (v_encounter_new_mag_sc^2) - (v_encounter_old_mag_remus^2))/(-2*v_encounter_new_mag_sc*v_encounter_old_mag_remus));
FPA_encounter_new_sc = kappa_angle - FPA_encounter_old_remus;
TA_encounter_new_sc = atand(((r_encounter*(v_encounter_new_mag_sc^2)/Gm_jupiter)*cosd(FPA_encounter_new_sc)*sind(FPA_encounter_new_sc))/(((r_encounter*(v_encounter_new_mag_sc^2)/Gm_jupiter)*(cosd(FPA_encounter_new_sc))^2)-1));

fprintf("delta_v_eq_mag: %.4e km/s\n", delta_v_eq_mag);
fprintf("nu_angle: %.4e deg\n", nu_angle);
fprintf("rho_angle: %.4e deg\n", rho_angle);
fprintf("rho_angle_check: %.4e deg\n", rho_angle_check);
fprintf("sigma_angle: %.4e deg\n", sigma_angle);
fprintf("v_encounter_new_mag_sc: %.4e km/s\n", v_encounter_new_mag_sc);
fprintf("kappa_angle: %.4e deg\n", kappa_angle);
fprintf("kappa_angle_check: %.4e deg\n", kappa_angle_check);
fprintf("FPA_encounter_new_sc: %.4e deg\n", FPA_encounter_new_sc);
fprintf("TA_encounter_new_sc: %.4e deg\n", TA_encounter_new_sc);




% plot function
function plot_eph(e,p,AOP, start_range, end_range)
    %orbit distance
    ta = [start_range:0.1:end_range];
    rmag = p ./ (1+e*cosd(ta-AOP));
    rplot(1,:) = rmag.*cosd(ta);
    rplot(2,:) = rmag.*sind(ta);
    
    plot(0,0,'.','MarkerSize',10,'Color','k','HandleVisibility','off'), hold on %Primary Body
    [~,indx] = min(abs(ta-AOP));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    [~,indx] = min(abs(ta-(180+AOP)));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit

    title("Lunar Free-Return Maneuver (Moon-centered View)-Lillian Shido")
    xlabel("Distance [km]")
    ylabel("Distance [km]")
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end