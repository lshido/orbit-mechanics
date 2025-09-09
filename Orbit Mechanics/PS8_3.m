a_earth = 149597898;
a_venus = 108207284;
Gm_earth = 398600.4415;
Gm_venus = 324858.59882646;
R_earth = 6378.1363; % [km]
R_venus = 6051.9;
Gm_sun = 132712440017.99;



% Problem 3a: Calculate transfer orbit
fprintf("----Problem 3a: Transfer Ellipse-----------\n");
rp_transfer = a_venus;
ra_transfer = a_earth;
a_transfer = (rp_transfer+ra_transfer)/2;
e_transfer = 1 - (rp_transfer/a_transfer);
vp_transfer = sqrt((2*(Gm_sun)/rp_transfer)-((Gm_sun)/a_transfer));
va_transfer = sqrt((2*(Gm_sun)/ra_transfer)-((Gm_sun)/a_transfer));
period_transfer = 2*pi*sqrt((a_transfer^3)/Gm_sun);
TOF_earth_venus = 1/2*period_transfer;
fprintf("rp_transfer: %.4e km\n", rp_transfer);
fprintf("ra_transfer: %.4e km\n", ra_transfer);
fprintf("a_transfer: %.4e km\n", a_transfer);
fprintf("e_transfer: %.4e km\n", e_transfer);
fprintf("vp_transfer: %.4e km/s\n", vp_transfer);
fprintf("va_transfer: %.4e km/s\n", va_transfer);
fprintf("period_transfer: %.4e sec\n", period_transfer);
fprintf("period_transfer: %.4e days\n", period_transfer/3600/24);
fprintf("period_transfer: %.4e years\n", period_transfer/3600/24/365);
fprintf("TOF_earth_venus: %.4e sec\n", TOF_earth_venus);
fprintf("TOF_earth_venus: %.4e days\n", TOF_earth_venus/3600/24);
fprintf("TOF_earth_venus: %.4e years\n", TOF_earth_venus/3600/24/365);

% Problem 3a: Calculate phase angle
fprintf("----Problem 3a: Phase Angle-----------\n");
mean_motion = sqrt((Gm_sun+Gm_venus)/a_venus^3);
phase_angle = 180 - rad2deg(mean_motion*TOF_earth_venus);
fprintf("mean_motion: %.4e km/s\n", mean_motion);
fprintf("phase_angle: %.4e deg\n", phase_angle);

% Problem 3a: Calculate planetary velocities
fprintf("----Problem 3a: Planetary Velocities-----------\n");
v_venus = sqrt((Gm_sun+Gm_venus)/rp_transfer);
v_earth = sqrt((Gm_sun+Gm_earth)/ra_transfer);
fprintf("v_venus: %.4e km/s\n", v_venus);
fprintf("v_earth: %.4e km/s\n", v_earth);

% Problem 3b: Calculate Departure from Earth
fprintf("----Problem 3b: Calculate Departure from Earth-----------\n");
r_parking_earth = R_earth + 210;
v_parking_earth = sqrt(Gm_earth/r_parking_earth);
v_inf_dep = va_transfer - v_earth;
delta_v_dep = sqrt((v_inf_dep^2) + (2*(Gm_earth)/r_parking_earth)) - v_parking_earth;
fprintf("v_parking_earth: %.4e km/s\n", v_parking_earth);
fprintf("r_parking_earth: %.4e km\n", r_parking_earth);
fprintf("v_inf_dep: %.4e km/s\n", v_inf_dep);
fprintf("delta_v_dep: %.4e km/s\n", delta_v_dep);

% Problem 3c: Calculate Arrival Conditions at Venus
fprintf("----Problem 3c: Calculate Arrival Conditions at Venus-----------\n");
r_old = a_venus;
rp_old = a_venus;
ra_old = a_earth;
v_old = vp_transfer;
energy_old = ((v_old^2)/2)-(Gm_sun/r_old);
energy_old_check = -Gm_sun/(2*a_transfer);
FPA_old = 0;
TA_old = 0;
fprintf("r_old: %.4e km\n", r_old);
fprintf("rp_old: %.4e km\n", rp_old);
fprintf("ra_old: %.4e km\n", ra_old);
fprintf("v_old: %.4e km/s\n", v_old);
fprintf("energy_old: %.4e km^2/s^2\n", energy_old);
fprintf("energy_old_check: %.4e km^2/s^2\n", energy_old_check);
fprintf("FPA_old: %.4e deg\n", FPA_old);
fprintf("TA_old: %.4e deg\n", TA_old);

% Problem 3c: Calculate orbit insertion at Venus
fprintf("----Problem 3c: Calculate orbit insertion at Venus-----------\n");
r_parking_venus = R_venus + 2000;
v_parking_venus = sqrt(Gm_venus/r_parking_venus);
energy_parking_venus = ((v_parking_venus^2)/2)-(Gm_venus/r_parking_venus);
v_inf_arr = vp_transfer - v_venus;
delta_v_arr = sqrt((v_inf_arr^2) + (2*(Gm_venus)/r_parking_venus)) - v_parking_venus;
fprintf("r_parking_venus: %.4e km\n", r_parking_venus);
fprintf("v_parking_venus: %.4e km/s\n", v_parking_venus);
fprintf("energy_parking_venus: %.4e km^2/s^2\n", energy_parking_venus);
fprintf("v_inf_arr: %.4e km/s\n", v_inf_arr);
fprintf("delta_v_arr: %.4e km/s\n", delta_v_arr);

% Problem 3c: Total Delta V
fprintf("----Problem 3c: Total Delta V-----------\n");
delta_v_total = delta_v_dep + delta_v_arr;
fprintf("delta_v_total 1: %.4e km/s\n", delta_v_arr);
fprintf("delta_v_total 2: %.4e km/s\n", delta_v_dep);
fprintf("delta_v_total: %.4e km/s\n", delta_v_total);

% Problem 3d: Flyby Hyperbola Characteristics
fprintf("----Problem 3d: Flyby Hyperbola Characteristics-----------\n");
energy_hyp_flyby = v_inf_arr^2/2;
a_hyp_flyby = Gm_venus/(2*energy_hyp_flyby);
e_hyp_flyby = (r_parking_venus/a_hyp_flyby) + 1;
turn_angle = asind(1/e_hyp_flyby);
FBA = 2*asind(1/e_hyp_flyby);
v_inf_arr_plus_vector = [ v_inf_arr*sind(FBA) v_inf_arr*cosd(FBA)];
kappa_angle = (180-FBA)/2;
alpha_angle = 180-kappa_angle;
delta_v_eq = 2*v_inf_arr*cosd(kappa_angle);
fprintf("energy_hyp_flyby: %.4e km^2/s^2\n", energy_hyp_flyby);
fprintf("a_hyp_flyby: %.4e km\n", a_hyp_flyby);
fprintf("e_hyp_flyby: %.4e\n", e_hyp_flyby);
fprintf("turn_angle: %.4e deg\n", turn_angle);
fprintf("FBA: %.4e deg\n", FBA);
fprintf("v_inf_arr_plus_vector: %.4e km/s r^ %.4e km/s theta^\n", v_inf_arr_plus_vector);
fprintf("kappa_angle: %.4e deg\n", kappa_angle);
fprintf("alpha_angle: %.4e deg\n", alpha_angle);
fprintf("delta_v_eq: %.4e km/s\n", delta_v_eq);

% Problem 3d: New Heliocentric Characteristics:
fprintf("----Problem 3d: New Heliocentric Characteristics-----------\n");
v_new = sqrt((v_old^2)+(delta_v_eq^2)-(2*v_old*delta_v_eq*cosd(kappa_angle)));
FPA_new = asind((delta_v_eq*sind(kappa_angle))/v_new);
TA_new = atand(((rp_transfer*(v_new^2)/Gm_sun)*cosd(FPA_new)*sind(FPA_new))/(((rp_transfer*(v_new^2)/Gm_sun)*(cosd(FPA_new))^2)-1));
e_new = sqrt(((((rp_transfer*v_new^2)/Gm_sun)-1)^2)*((cosd(FPA_new))^2)+(sind(FPA_new))^2);
energy_new = ((v_new^2)/2)-(Gm_sun/rp_transfer);
a_new = -Gm_sun/(2*energy_new);
rp_new = a_new*(1 - e_new);
ra_new = a_new*(1 + e_new);
delta_small_omega = -TA_new + TA_old;
p_new = a_new*(1-(e_new)^2);

fprintf("v_new: %.4e km/s\n", v_new);
fprintf("FPA_new: %.4e km/s\n", FPA_new);
fprintf("TA_new: %.4e km/s\n", TA_new);
fprintf("e_new: %.4e\n", e_new);
fprintf("energy_new: %.4e km/s\n", energy_new);
fprintf("a_new: %.4e km\n", a_new);
fprintf("rp_new: %.4e km\n", rp_new);
fprintf("ra_new: %.4e km\n", ra_new);
fprintf("delta_small_omega: %.4e deg\n", delta_small_omega);




% Earth orbit
plot_eph(0, a_earth, 0, 360);

% Venus orbit
plot_eph(0, a_venus, 0, 360);

% Mercury orbit
a_mercury = 57909101;
plot_eph(0, a_mercury, 0, 360);

% heliocentric orbit
plot_eph(e_new, p_new, delta_small_omega, 360)

% transfer orbit
p_transfer = a_transfer*(1-(e_transfer)^2);
plot_eph(e_transfer, p_transfer, 0, 360)

% plot function
function plot_eph(e,p,AOP, periodic_range, orbit_name)
    %orbit distance
    ta = [0:0.1:periodic_range];
    rmag = p ./ (1+e*cosd(ta-AOP));
    rplot(1,:) = rmag.*cosd(ta);
    rplot(2,:) = rmag.*sind(ta);
    
    plot(0,0,'.','MarkerSize',25,'Color','k','HandleVisibility','off'), hold on %Primary Body
    [~,indx] = min(abs(ta-AOP));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    [~,indx] = min(abs(ta-(180+AOP)));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit

    title("Transfer Plan of Spacecraft to Venus-Lillian Shido")
    xlabel("Distance [km]")
    ylabel("Distance [km]")
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end

% plot function
function plot_eph_transfer(e,p,AOP)
    %orbit distance
    ta = [0:0.1:180];
    rmag = p ./ (1+e*cosd(ta-AOP));
    rplot(1,:) = rmag.*cosd(ta);
    rplot(2,:) = rmag.*sind(ta);
    
    plot(0,0,'.','MarkerSize',25,'Color','r','HandleVisibility','off'), hold on %Primary Body
    [~,indx] = min(abs(ta-AOP));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    [~,indx] = min(abs(ta-(180+AOP)));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end