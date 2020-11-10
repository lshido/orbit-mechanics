Gm_earth = 398600.4415;
R_earth = 6378.1363; % [km]
Gm_moon = 4902.8005821478;
a_earth_moon = +384400;
R_moon = 1738.2;


% Problem 1a: Calculate dark-side arrival
fprintf("-----------Problem 1a: Calculate dark-side arrival---------\n")
rp_transfer = R_earth + 190;
ra_transfer = 3.8440e5; % km
va_transfer = 0.18666; % km/s
v_inf_arr = 0.83789; % km/s
FBA = 103.03; % deg
fprintf("rp_transfer: %.4e km\n", rp_transfer);
fprintf("ra_transfer: %.4e km\n", ra_transfer);
fprintf("FBA: %.4e deg\n", FBA);

% Problem 1a: Orbital Characteristics after Flyby
fprintf("-----------Problem 1a: Orbital Characteristics after Flyby---------\n")
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

% % iterative solver for a_free_return - Lillian Shido
% fprintf("-----------Problem 1a: Iteratively solve for a_transfer---------\n")
% r = ra_transfer;
% TA = 173.8;
% r_error = 1;
% a_guess = ra_transfer+1e5;
% stepsize = 0.00001e5;

% while abs(r_error) > 0.01
%     if r_error > 0
%         a_guess = a_guess - stepsize;
%     else
%         a_guess = a_guess + stepsize;
%     end
%     r_calced = ((rp_transfer/a_guess)*((2*a_guess)-rp_transfer))/(1+((1-(rp_transfer/a_guess))*cosd(TA)));
%     r_error = r_calced - r;
% end
% fprintf("Calculated a_free_return: %.4e km\n", a_guess);

% Problem 1b: Free-return Transfer Characteristics
fprintf("-----------Problem 1b: Free-return Transfer Characteristics---------\n")
r_free_return = a_earth_moon;
TA = 173.8;
a_free_return = 2.3450e5;
e_free_return = 1 - rp_transfer/a_free_return;
p_free_return = a_free_return*(1-e_free_return^2);
period_transfer = 2*pi*sqrt((a_free_return^3)/Gm_earth);
rp_free_return = a_free_return*(1-e_free_return);
ra_free_return = a_free_return*(1+e_free_return);
energy_free_return = -Gm_earth/(2*a_free_return);
EA_free_return = 2*atan(sqrt((1-e_free_return)/(1+e_free_return))*tand(TA/2));
TOF_outbound = sqrt((a_free_return^3)/Gm_earth)*(EA_free_return - e_free_return*sin(EA_free_return));
phi_angle = TA - rad2deg(sqrt(Gm_earth/(a_free_return^3))*TOF_outbound);
fprintf("a_free_return: %.4e\n", a_free_return);
fprintf("p_free_return: %.4e\n", p_free_return);
fprintf("e_free_return: %.4e\n", e_free_return);
fprintf("rp_free_return: %.4e km\n", rp_free_return);
fprintf("ra_free_return: %.4e km\n", ra_free_return);
fprintf("energy_free_return: %.4e km^2/s^2\n", energy_free_return);
fprintf("period_transfer: %.4e sec\n", period_transfer);
fprintf("period_transfer: %.4e days\n", period_transfer/3600/24);
fprintf("period_transfer: %.4e years\n", period_transfer/3600/24/365);
fprintf("EA_free_return: %.4e rad\n", EA_free_return);
fprintf("EA_free_return: %.4e deg\n", rad2deg(EA_free_return));
fprintf("TOF_outbound: %.4e sec\n", TOF_outbound);
fprintf("TOF_outbound: %.4e hours\n", TOF_outbound/3600);
fprintf("TOF_outbound: %.4e days\n", TOF_outbound/3600/24);
fprintf("phi_angle: %.4e deg\n", phi_angle);

% Problem 1c: Free-return Lunar Arrival
fprintf("-----------Problem 1c: Free-return Lunar Arrival---------\n")
v_free_return_old_mag = sqrt(2*(energy_free_return + (Gm_earth/r_free_return)));
FPA_free_return_old = acosd(sqrt(p_free_return*Gm_earth)/(r_free_return*v_free_return_old_mag));
v_free_return_old_rhat = sqrt((v_free_return_old_mag^2)-(p_free_return*Gm_earth/r_free_return^2));
v_free_return_old_thetahat = sqrt(p_free_return*Gm_earth)/r_free_return;
v_free_return_old_rhat_check = v_free_return_old_mag*sind(FPA_free_return_old);
v_free_return_old_thetahat_check = v_free_return_old_mag*cosd(FPA_free_return_old);
v_moon_mag = sqrt((Gm_earth+Gm_moon)/a_earth_moon);
v_free_return_old = [ v_free_return_old_rhat v_free_return_old_thetahat ];
v_moon = [ 0 v_moon_mag ];
v_inf_arr_free_return = v_free_return_old - v_moon;
v_inf_arr_free_return_mag = norm(v_inf_arr_free_return);
turn_angle = asind(v_free_return_old_mag*sind(FPA_free_return_old)/v_inf_arr_free_return_mag);
FBA_free_return_hyp = 2*turn_angle;
e_free_return_hyp = 1/(sind(FBA_free_return_hyp)/2);
energy_free_return_hyp = (v_inf_arr_free_return_mag^2)/2;
a_free_return_hyp = Gm_moon/(2*energy_free_return_hyp);
p_free_return_hyp = a_free_return_hyp*((e_free_return_hyp^2)-1);
r_pass_free_return = a_free_return_hyp*(e_free_return_hyp-1);
passing_altitude = r_pass_free_return - R_moon;
fprintf("v_free_return_old_mag: %.4e km/s\n", v_free_return_old_mag);
fprintf("v_free_return_old_rhat: %.4e km/s\n", v_free_return_old_rhat);
fprintf("v_free_return_old_thetahat: %.4e km/s\n", v_free_return_old_thetahat);
fprintf("FPA_free_return_old: %.4e deg\n", FPA_free_return_old);
fprintf("v_free_return_old_rhat_check: %.4e km/s\n", v_free_return_old_rhat_check);
fprintf("v_free_return_old_thetahat_check: %.4e km/s\n", v_free_return_old_thetahat_check);
fprintf("v_moon_mag: %.4e km/s\n", v_moon_mag);
fprintf("v_inf_arr_free_return: %.4e km/s r^ %.4e km/s theta^\n", v_inf_arr_free_return);
fprintf("v_inf_arr_free_return_mag: %.4e km/s\n", v_inf_arr_free_return_mag);
fprintf("turn_angle: %.4e deg\n", turn_angle);
fprintf("FBA_free_return_hyp: %.4e deg\n", FBA_free_return_hyp);
fprintf("e_free_return_hyp: %.4e\n", e_free_return_hyp);
fprintf("energy_free_return_hyp: %.4e km^2/s^2\n", energy_free_return_hyp);
fprintf("a_free_return_hyp: %.4e km\n", a_free_return_hyp);
fprintf("r_pass_free_return: %.4e km\n", r_pass_free_return);
fprintf("passing_altitude: %.4e km\n", passing_altitude);


fprintf("-----------Problem 1c: Free-return Lunar Departurel--------\n")
delta_v_eq_free_return_mag = 2*v_inf_arr_free_return_mag*sind(FBA_free_return_hyp/2);
delta_v_eq_free_return = [ -delta_v_eq_free_return_mag 0 ];
v_free_return_new = v_free_return_old + delta_v_eq_free_return;
v_free_return_new_mag = norm(v_free_return_new);
FPA_free_return_new = atand(v_free_return_new(1)/v_free_return_new(2));
TA_free_return_new = 180 + atand(((r_free_return*(v_free_return_new_mag^2)/Gm_earth)*cosd(FPA_free_return_new)*sind(FPA_free_return_new))/(((r_free_return*(v_free_return_new_mag^2)/Gm_earth)*(cosd(FPA_free_return_new))^2)-1));
energy_free_return_new = ((v_free_return_new_mag^2)/2)-(Gm_earth/r_free_return);
a_free_return_new = -Gm_earth/(2*energy_free_return_new);
e_free_return_new = sqrt(((((r_free_return*v_free_return_new_mag^2)/Gm_earth)-1)^2)*((cosd(FPA_free_return_new))^2)+(sind(FPA_free_return_new))^2);
period_free_return_new = 2*pi*sqrt((a_free_return_new^3)/Gm_earth);
rp_free_return_new = a_free_return_new*(1-e_free_return_new);
ra_free_return_new = a_free_return_new*(1+e_free_return_new);
delta_small_omega_free_return_new = -TA_free_return_new + TA;

fprintf("delta_v_eq_free_return_mag: %.4e km/s\n", delta_v_eq_free_return_mag);
fprintf("delta_v_eq_free_return: %.4e km/s r^ %.4e km/s theta^\n", delta_v_eq_free_return);
fprintf("v_free_return_new: %.4e km/s r^ %.4e km/s theta^\n", v_free_return_new);
fprintf("v_free_return_new_mag: %.4e km/s\n", v_free_return_new_mag);
fprintf("FPA_free_return_new: %.4e deg\n", FPA_free_return_new);
fprintf("TA_free_return_new: %.4e deg\n", TA_free_return_new);
fprintf("energy_free_return_new: %.4e km^2/s^2\n", energy_free_return_new);
fprintf("a_free_return_new: %.4e km\n", a_free_return_new);
fprintf("e_free_return_new: %.4e km\n", e_free_return_new);
fprintf("period_free_return_new: %.4e sec\n", period_free_return_new);
fprintf("period_free_return_new: %.4e days\n", period_free_return_new/3600/24);
fprintf("period_free_return_new: %.4e years\n", period_free_return_new/3600/24/365);
fprintf("rp_free_return_new: %.4e km\n", rp_free_return_new);
fprintf("ra_free_return_new: %.4e km\n", ra_free_return_new);
fprintf("delta_small_omega_free_return_new: %.4e deg\n", delta_small_omega_free_return_new);


fprintf("-----------Problem 1c: Convert to VNB and calc alpha--------\n")
TM_vnb_to_polar = [ sind(FPA_free_return_old) cosd(FPA_free_return_old); cosd(FPA_free_return_old) -sind(FPA_free_return_old)];
TM_polar_to_VNB = transpose(TM_vnb_to_polar);
delta_v_eq_free_return_vnb = delta_v_eq_free_return*TM_polar_to_VNB;
kappa_angle_new = 180 - (90+FPA_free_return_old);
alpha_angle_new = 180-kappa_angle_new;
fprintf("delta_v_eq_free_return_vnb: %.4e km/s V^ %.4e km/s B^\n", delta_v_eq_free_return_vnb);
fprintf("kappa_angle_new: %.4e deg\n", kappa_angle_new);
fprintf("alpha_angle_new: %.4e deg\n", alpha_angle_new);


% % Earth-centered view
% % parking orbit
% plot_eph(0, rp_transfer, 0, 0, 360);
% % outbound
% plot_eph(e_free_return, p_free_return, 6.2, 0, 180)
% % return
% plot_eph(e_free_return, p_free_return, -6.2, 180, 360)
% % moon orbit
% plot_eph(0, 384400,0, 0, 360)

% luna-centered view
plot_eph(e_free_return_hyp, p_free_return_hyp, 0, 0, 360)


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