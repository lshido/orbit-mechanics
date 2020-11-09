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

% Problem 1a: Free-return Transfer Characteristics
fprintf("-----------Problem 1a: Free-return Transfer Characteristics---------\n")
a_free_return = 2.3450e5;
e_free_return = 1 - rp_transfer/a_guess;
period_transfer = 2*pi*sqrt((a_free_return^3)/Gm_earth);
rp_free_return = a_free_return*(1-e_free_return);
ra_free_return = a_free_return*(1+e_free_return);
energy_free_return = -Gm_earth/(2*a_free_return);
EA_free_return = 2*atan(sqrt((1-e_free_return)/(1+e_free_return))*tand(TA/2));
TOF_outbound = sqrt((a_free_return^3)/Gm_earth)*(EA_free_return - e_free_return*sin(EA_free_return));
phi_angle = TA - rad2deg(sqrt(Gm_earth/(a_free_return^3))*TOF_outbound);
fprintf("a_free_return: %.4e\n", a_free_return);
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


p=1.9548e5*(1-0.96640^2);
plot_eph(0.96640, p, 0);
p_free_return = a_free_return*(1-e_free_return^2);
plot_eph(e_free_return, p_free_return, 0)
plot_eph(0, 384400,0)

% plot function
function plot_eph(e,p,AOP)
    %orbit distance
    ta = [0:0.1:360];
    rmag = p ./ (1+e*cosd(ta-AOP));
    rplot(1,:) = rmag.*cosd(ta);
    rplot(2,:) = rmag.*sind(ta);
    
    plot(0,0,'.','MarkerSize',25,'Color','k','HandleVisibility','off'), hold on %Primary Body
    [~,indx] = min(abs(ta-AOP));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    [~,indx] = min(abs(ta-(180+AOP)));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit

    title("Flyby of Spacecraft at Jupiter-Lillian Shido")
    xlabel("Distance [km]")
    ylabel("Distance [km]")
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end