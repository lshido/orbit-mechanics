Gm_sun = 132712440017.99;
Gm_venus = 324858.59882646;
R_venus = 6051.9;

earth_launch = 2453221.0;
first_venus_flyby = 2454033.0;
first_mercury_flyby = 2454480.0;
mercury_orbit_insertion = 2455639.0;

fprintf("-----------Problem 2a: TOFs--------\n")
earth_launch_to_first_venus_flyby = first_venus_flyby - earth_launch;
earth_launch_to_first_mercury_flyby = first_mercury_flyby - earth_launch;
earth_launch_to_mercury_orbit_insertion = mercury_orbit_insertion - earth_launch;
fprintf("earth_launch_to_first_venus_flyby: %.4e days\n", earth_launch_to_first_venus_flyby);
fprintf("earth_launch_to_first_venus_flyby: %.4e years\n", earth_launch_to_first_venus_flyby/365);
fprintf("earth_launch_to_first_mercury_flyby: %.4e days\n", earth_launch_to_first_mercury_flyby);
fprintf("earth_launch_to_first_mercury_flyby: %.4e years\n", earth_launch_to_first_mercury_flyby/365);
fprintf("earth_launch_to_mercury_orbit_insertion: %.4e days\n", earth_launch_to_mercury_orbit_insertion);
fprintf("earth_launch_to_mercury_orbit_insertion: %.4e years\n", earth_launch_to_mercury_orbit_insertion/365);

fprintf("-----------Problem 2b: Compute TA at venus encounter--------\n")
rp_Hohmann = 7.4800e7;
ra_Hohmann = 1.4960e8;
r_encounter = 1.0821e8;
venus_au = r_encounter/ra_Hohmann;
a_Hohmann = (rp_Hohmann + ra_Hohmann)/2;
e_Hohmann = 1 - (rp_Hohmann/a_Hohmann);
p_Hohmann = a_Hohmann*(1-(e_Hohmann^2));
energy_Hohmann = -Gm_sun/(2*a_Hohmann);
period_Hohmann = 2*pi*sqrt((a_Hohmann^3)/Gm_sun);
TA_encounter = -acosd((p_Hohmann/(r_encounter*e_Hohmann))-(1/e_Hohmann));
EA_encounter = 2*atand(sqrt((1-e_Hohmann)/(1+e_Hohmann))*tand(TA_encounter/2));
time_to_venus = sqrt((a_Hohmann^3)/Gm_sun)*(deg2rad(EA_encounter) - (e_Hohmann*sind(EA_encounter)));
time_to_earth = period_Hohmann/2;
TOF_earth_to_venus = (period_Hohmann + time_to_venus) - time_to_earth;
fprintf("venus_au: %.4e AU\n", venus_au);
fprintf("rp_Hohmann: %.4e km\n", rp_Hohmann);
fprintf("ra_Hohmann: %.4e km\n", ra_Hohmann);
fprintf("a_Hohmann: %.4e km\n", a_Hohmann);
fprintf("e_Hohmann: %.4e\n", e_Hohmann);
fprintf("p_Hohmann: %.4e\n", p_Hohmann);
fprintf("energy_Hohmann: %.4e km^2/s^2\n", energy_Hohmann);
fprintf("TA_encounter: %.4e deg\n", TA_encounter);
fprintf("EA_encounter: %.4e deg\n", EA_encounter);
fprintf("EA_encounter: %.4e rad\n", deg2rad(EA_encounter));
fprintf("time_to_venus: %.4e sec\n", time_to_venus);
fprintf("time_to_venus: %.4e days\n", time_to_venus/3600/24);
fprintf("time_to_earth: %.4e sec\n", time_to_earth);
fprintf("TOF_earth_to_venus: %.4e sec\n", TOF_earth_to_venus);
fprintf("TOF_earth_to_venus: %.4e days\n", TOF_earth_to_venus/3600/24);
fprintf("period_Hohmann: %.4e sec\n", period_Hohmann);
fprintf("period_Hohmann: %.4e days\n", period_Hohmann/3600/24);
fprintf("period_Hohmann: %.4e years\n", period_Hohmann/3600/24/365);

fprintf("-----------Problem 2c: Compute post-Venus Encounter--------\n")
rp_venus = R_venus + 2990;
v_venus_mag = sqrt((Gm_sun + Gm_venus)/r_encounter);
v_venus = [ 0 v_venus_mag ];
v_encounter_old_mag = sqrt(2*(energy_Hohmann+(Gm_sun/r_encounter)));
FPA_encounter_old = -acosd(sqrt(p_Hohmann*Gm_sun)/(r_encounter*v_encounter_old_mag));
v_encounter_old = [ v_encounter_old_mag*sind(FPA_encounter_old) v_encounter_old_mag*cosd(FPA_encounter_old) ];
v_encounter_inf = v_encounter_old - v_venus;
v_encounter_inf_mag = norm(v_encounter_inf);
nu_encounter_hyp = acosd(((v_encounter_old_mag^2) - (v_venus_mag^2) - (v_encounter_inf_mag^2))/(-2*v_venus_mag*v_encounter_inf_mag));
a_encounter_hyp = Gm_venus/(v_encounter_inf_mag^2);
e_encounter_hyp = (rp_venus/a_encounter_hyp)+1;
turn_angle_encounter = asind(1/e_encounter_hyp);
FBA_encounter_hyp = 2*turn_angle_encounter;
delta_v_encounter_mag = 2*v_encounter_inf_mag*sind(FBA_encounter_hyp/2);
rho_angle = nu_encounter_hyp - FBA_encounter_hyp;
v_encounter_new_mag = sqrt((v_venus_mag^2)+(v_encounter_inf_mag^2)-(2*v_venus_mag*v_encounter_inf_mag*cosd(rho_angle)));
FPA_encounter_new = -asind((v_encounter_inf_mag*sind(rho_angle))/v_encounter_new_mag);
v_encounter_new = [ v_encounter_new_mag*sind(FPA_encounter_new) v_encounter_new_mag*cosd(FPA_encounter_new) ];
TA_encounter_new = 180 + atand(((r_encounter*(v_encounter_new_mag^2)/Gm_sun)*cosd(FPA_encounter_new)*sind(FPA_encounter_new))/(((r_encounter*(v_encounter_new_mag^2)/Gm_sun)*(cosd(FPA_encounter_new))^2)-1));
e_encounter_new = sqrt(((((r_encounter*v_encounter_new_mag^2)/Gm_sun)-1)^2)*((cosd(FPA_encounter_new))^2)+(sind(FPA_encounter_new))^2);
energy_encounter_new = ((v_encounter_new_mag^2)/2)-(Gm_sun/r_encounter);
a_encounter_new = -Gm_sun/(2*energy_encounter_new);
rp_encounter_new = a_encounter_new*(1-e_encounter_new);
ra_encounter_new = a_encounter_new*(1+e_encounter_new);
period_encounter_new = 2*pi*sqrt((a_encounter_new^3)/Gm_sun);
delta_small_omega_encounter = -TA_encounter_new + TA_encounter;
delta_v_encounter = v_encounter_new - v_encounter_old;
beta_angle = acosd(((v_encounter_new_mag^2) - (v_encounter_old_mag^2) - (delta_v_encounter_mag^2))/(-2*v_encounter_old_mag*delta_v_encounter_mag));
delta_FPA_encounter = FPA_encounter_new - FPA_encounter_old;
beta_angle_check = asind((v_encounter_new_mag*sind(delta_FPA_encounter))/delta_v_encounter_mag);
alpha_angle_encounter = -(180 - beta_angle);
v_encounter_new_vnb = [ v_encounter_new_mag*cosd(delta_FPA_encounter) v_encounter_new_mag*sind(delta_FPA_encounter) ];
delta_v_encounter_vnb = [ delta_v_encounter_mag*cosd(alpha_angle_encounter) delta_v_encounter_mag*sind(alpha_angle_encounter) ];

fprintf("rp_venus: %.4e km\n", rp_venus);
fprintf("rp_venus radii: %.4e km\n", rp_venus/R_venus);
fprintf("v_venus_mag: %.4e km/s\n", v_venus_mag);
fprintf("v_encounter_old_mag: %.4e km/s\n", v_encounter_old_mag);
fprintf("FPA_encounter_old: %.4e deg\n", FPA_encounter_old);
fprintf("v_venus: %.4e km/s r^ %.4e km/s theta^\n", v_venus);
fprintf("v_encounter_old: %.4e km/s r^ %.4e km/s theta^\n", v_encounter_old);
fprintf("v_encounter_inf: %.4e km/s r^ %.4e km/s theta^\n", v_encounter_inf);
fprintf("v_encounter_inf_mag: %.4e km/s\n", v_encounter_inf_mag);
fprintf("nu_encounter_hyp: %.4e deg\n", nu_encounter_hyp);
fprintf("delta_v_encounter_mag: %.4e km/s\n", delta_v_encounter_mag);
fprintf("a_encounter_hyp: %.4e km\n", a_encounter_hyp);
fprintf("e_encounter_hyp: %.4e\n", e_encounter_hyp);
fprintf("turn_angle_encounter: %.4e deg\n", turn_angle_encounter);
fprintf("FBA_encounter_hyp: %.4e deg\n", FBA_encounter_hyp);
fprintf("rho_angle: %.4e deg\n", rho_angle);
fprintf("v_encounter_new_mag: %.4e km/s\n", v_encounter_new_mag);
fprintf("FPA_encounter_new: %.4e deg\n", FPA_encounter_new);
fprintf("v_encounter_new: %.4e km/s r^ %.4e km/s theta^\n", v_encounter_new);
fprintf("TA_encounter_new: %.4e deg\n", TA_encounter_new);
fprintf("e_encounter_new: %.4e\n", e_encounter_new);
fprintf("energy_encounter_new: %.4e km^2/s^2\n", energy_encounter_new);
fprintf("a_encounter_new: %.4e km\n", a_encounter_new);
fprintf("rp_encounter_new: %.4e km\n", rp_encounter_new);
fprintf("ra_encounter_new: %.4e km\n", ra_encounter_new);
fprintf("period_encounter_new: %.4e sec\n", period_encounter_new);
fprintf("period_encounter_new: %.4e days\n", period_encounter_new/3600/24);
fprintf("period_encounter_new: %.4e years\n", period_encounter_new/3600/24/365);
fprintf("delta_small_omega_encounter: %.4e deg\n", delta_small_omega_encounter);
fprintf("delta_v_encounter: %.4e km/s r^ %.4e km/s theta^\n", delta_v_encounter);
fprintf("beta_angle: %.4e deg\n", beta_angle);
fprintf("delta_FPA_encounter: %.4e deg\n", delta_FPA_encounter);
fprintf("beta_angle_check: %.4e deg\n", beta_angle_check);
fprintf("alpha_angle_encounter: %.4e deg\n", alpha_angle_encounter);
fprintf("v_encounter_new_vnb: %.4e km/s V^ %.4e km/s B^\n", v_encounter_new_vnb);
fprintf("delta_v_encounter_vnb: %.4e km/s V^ %.4e km/s B^\n", delta_v_encounter_vnb);


% old orbit
plot_eph(e_Hohmann, p_Hohmann, 0, 0, 360)

% new orbit
p_encounter_new = a_encounter_new*(1 - e_encounter_new^2);
plot_eph(e_encounter_new, p_encounter_new, delta_small_omega_encounter, 0, 360)

% venus orbit
% plot_eph(0, r_encounter, 0, 0, 360)

% earth orbit
% plot_eph(0, ra_Hohmann, 0, 0, 360)

% mercury orbit
r_mercury = 57909101;
plot_eph(0, r_mercury, 0, 0, 360)

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
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','on') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit

    title("Non-tangential Venus Flyby (Heliocentric View)-Lillian Shido")
    xlabel("Distance [km]")
    ylabel("Distance [km]")
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end