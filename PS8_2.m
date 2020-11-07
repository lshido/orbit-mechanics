Gm_earth = 398600.4415;
R_earth = 6378.1363; % [km]
Gm_moon = 4902.8005821478;
a_earth_moon = +384400;
R_moon = 1738.2;


fprintf("----------Problem 2a-----------\n");
r_parking_earth = R_earth + 190;
r_parking_moon = R_moon + 200;
v_circular_earth = sqrt(Gm_earth/r_parking_earth);
v_circular_moon = sqrt(Gm_moon/r_parking_moon);
fprintf("r_parking_earth: %.4e km\n", r_parking_earth);
fprintf("r_parking_moon: %.4e km\n", r_parking_moon);
fprintf("v_circular_earth: %.4e km/s\n", v_circular_earth);
fprintf("v_circular_moon: %.4e km/s\n", v_circular_moon);

% calculate transfer ellipse
rp_transfer = r_parking_earth;
ra_transfer = a_earth_moon;
a_transfer = (rp_transfer+ra_transfer)/2;
e_transfer = 1 - (rp_transfer/a_transfer);
vp_transfer = sqrt((2*(Gm_earth)/rp_transfer)-((Gm_earth)/a_transfer));
va_transfer = sqrt((2*(Gm_earth)/ra_transfer)-((Gm_earth)/a_transfer));
fprintf("rp_transfer: %.4e km\n", rp_transfer);
fprintf("ra_transfer: %.4e km\n", ra_transfer);
fprintf("a_transfer: %.4e km\n", a_transfer);
fprintf("e_transfer: %.4e km\n", e_transfer);
fprintf("vp_transfer: %.4e km/s\n", vp_transfer);
fprintf("va_transfer: %.4e km/s\n", va_transfer);

% calculate delta_v_dep
v_inf_dep = vp_transfer;
delta_v_dep = sqrt((v_inf_dep^2) + (2*(Gm_earth)/r_parking_earth)) - v_circular_earth;
fprintf("v_inf_dep: %.4e km/s\n", v_inf_dep);
fprintf("delta_v_dep: %.4e km/s\n", delta_v_dep);

% calculate delta_v_arr
%v_moon_around_earth = sqrt((Gm_moon+Gm_earth)/a_earth_moon);
v_moon_around_earth = sqrt((Gm_earth)/a_earth_moon);
v_inf_arr = va_transfer - v_moon_around_earth;
delta_v_arr = sqrt((v_inf_arr^2) + (2*(Gm_moon)/r_parking_moon)) - v_circular_moon;
fprintf("v_moon_around_earth: %.4e km/s\n", v_moon_around_earth);
fprintf("v_inf_arr: %.4e km/s\n", v_inf_arr);
fprintf("delta_v_arr: %.4e km/s\n", delta_v_arr);

% calculate delta_v_total
delta_v_total = delta_v_arr + delta_v_dep;
fprintf("delta_v_total: %.4e km/s\n", delta_v_total);

% orbital parameters
fprintf("--Problem 2a: transfer ellipse orbital params-----------\n");
FPA_apoapsis = 0;
TA_apoapsis = 180;
period_transfer = 2*pi*sqrt((a_transfer^3)/Gm_earth);
TOF_earth_moon = 1/2*period_transfer;
energy_transfer = -Gm_earth/(2*a_transfer);
mean_motion = sqrt((Gm_earth+Gm_moon)/a_earth_moon^3);
phase_angle = 180 - rad2deg(mean_motion*TOF_earth_moon);
fprintf("ra_transfer: %.4e km\n", ra_transfer);
fprintf("rp_transfer: %.4e km\n", rp_transfer);
fprintf("va_transfer: %.4e km/s\n", va_transfer);
fprintf("FPA_apoapsis: %.4e deg\n", FPA_apoapsis);
fprintf("TA_apoapsis: %.4e deg\n", TA_apoapsis);
fprintf("a_transfer: %.4e km\n", a_transfer);
fprintf("e_transfer: %.4e\n", e_transfer);
fprintf("energy_transfer: %.4e km^2/s^2\n", energy_transfer);
fprintf("period_transfer: %.4e sec\n", period_transfer);
fprintf("period_transfer: %.4e days\n", period_transfer/3600/24);
fprintf("period_transfer: %.4e years\n", period_transfer/3600/24/365);
fprintf("TOF_earth_moon: %.4e sec\n", TOF_earth_moon);
fprintf("TOF_earth_moon: %.4e days\n", TOF_earth_moon/3600/24);
fprintf("TOF_earth_moon: %.4e years\n", TOF_earth_moon/3600/24/365);
fprintf("mean_motion: %.4e km/s\n", mean_motion);
fprintf("phase_angle: %.4e deg\n", phase_angle);


% Problem 2b: Characteristics of hyperbola for jupiter flyby
fprintf("--Problem 2b: Flyby Hyperbola Characteristics-----------\n");
energy_hyp_flyby = v_inf_arr^2/2;
a_hyp_flyby = Gm_moon/(2*energy_hyp_flyby);
e_hyp_flyby = (r_parking_moon/a_hyp_flyby) + 1;
turn_angle = asind(1/e_hyp_flyby);
FBA = 2*asind(1/e_hyp_flyby);
fprintf("energy_hyp_flyby: %.4e km^2/s^2\n", energy_hyp_flyby);
fprintf("a_hyp_flyby: %.4e km\n", a_hyp_flyby);
fprintf("e_hyp_flyby: %.4e\n", e_hyp_flyby);
fprintf("turn_angle: %.4e deg\n", turn_angle);
fprintf("FBA: %.4e deg\n", FBA);

% Problem 2b: Orbital Characteristics after Flyby
alpha_angle_flyby = (180-FBA)/2;
delta_v_flyby = 2*abs(v_inf_arr)*cosd(alpha_angle_flyby);
v_new_flyby = sqrt((va_transfer^2)+(delta_v_flyby^2)-2*va_transfer*delta_v_flyby*cosd(180-alpha_angle_flyby));
FPA_new_flyby = asind((delta_v_flyby*sind(180-alpha_angle_flyby))/v_new_flyby);
TA_new_flyby = atand(((ra_transfer*(v_new_flyby^2)/Gm_earth)*cosd(FPA_new_flyby)*sind(FPA_new_flyby))/(((ra_transfer*(v_new_flyby^2)/Gm_earth)*(cosd(FPA_new_flyby))^2)-1));
energy_new_flyby = ((v_new_flyby^2)/2)-(Gm_earth/ra_transfer);
a_new_flyby = Gm_earth/(2*energy_new_flyby);
v_new_flyby_rhat = v_new_flyby*sind(FPA_new_flyby);
v_new_flyby_thetahat = v_new_flyby*cosd(FPA_new_flyby);
e_flyby_new = sqrt(((((ra_transfer*v_new_flyby^2)/Gm_earth)-1)^2)*((cosd(FPA_new_flyby))^2)+(sind(FPA_new_flyby))^2);
rp_flyby_new = a_new_flyby*(e_flyby_new-1);
delta_small_omega = -TA_new_flyby + 180;
delta_v_flyby_vector = [ delta_v_flyby*cosd(alpha_angle_flyby) 0 delta_v_flyby*sind(alpha_angle_flyby) ];
fprintf("alpha_angle_flyby: %.4e deg\n", alpha_angle_flyby);
fprintf("delta_v_flyby: %.4e km/s\n", delta_v_flyby);
fprintf("v_new_flyby: %.4e km/s\n", v_new_flyby);
fprintf("FPA_new_flyby: %.4e km/s\n", FPA_new_flyby);
fprintf("TA_new_flyby: %.4e km/s\n", TA_new_flyby);
fprintf("energy_new_flyby: %.4e km/s\n", energy_new_flyby);
fprintf("a_new_flyby: %.4e km\n", a_new_flyby);
fprintf("v_new_flyby_rhat: %.4e km/s\n", v_new_flyby_rhat);
fprintf("v_new_flyby_thetahat: %.4e km/s\n", v_new_flyby_thetahat);
fprintf("e_flyby_new: %.4e\n", e_flyby_new);
fprintf("rp_flyby_new: %.4e km\n", rp_flyby_new);
fprintf("delta_small_omega: %.4e deg\n", delta_small_omega);
fprintf("delta_v_flyby_vector: %.4e km/s V^ %.4e km/s N^ %.4e km/s B^\n", delta_v_flyby_vector);

