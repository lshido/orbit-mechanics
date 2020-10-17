R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
a = 3*R_earth;
e=0.6;

r = 3*R_earth;
fprintf("radius at maneuver point: %.4e km\n", r);

%calc specific energy
specific_energy_old = -Gm_earth/(2*a);
fprintf("specific energy before maneuver: %.4e km^2/sec^2\n", specific_energy_old);

% calc v
v_old = sqrt(2*(specific_energy_old+(Gm_earth/r)));
fprintf("Velocity mag before maneuver: %.4e km/sec\n", v_old);

% calculate semilatus rectum
p_old = a*(1-(e^2));
fprintf("semilatus rectum: %.4e km\n", p_old);


% calc r_dot 
r_dot = -sqrt((v_old^2)-(p_old*Gm_earth/(r^2))); % negative because descending
fprintf("r_dot mag before maneuver: %.4e km/sec\n", r_dot);

% calc r_theta_dot
r_theta_dot = sqrt(p_old*Gm_earth)/r;
fprintf("r_theta_dot mag before maneuver: %.4e km/sec\n", r_theta_dot);

% calc FPA
FPA_old = atand(r_dot/r_theta_dot);
fprintf("FPA before maneuver: %.4e deg\n", FPA_old);

% calc old TA
TA_old = -acosd((p_old/(r*e))-(1/e)); % negative because descending
fprintf("True Anomaly before maneuver: %.4e deg\n", TA_old);
fprintf("True Anomaly before maneuver: %.4e rad\n", deg2rad(TA_old));

% express r and v in inertial frame
TM_inertial = [ cosd(TA_old) sind(TA_old); -sind(TA_old) cosd(TA_old)];
r_inertial = [r 0] * TM_inertial;
v_inertial_old = [r_dot r_theta_dot] * TM_inertial;
fprintf("r_inertial: %.4e km e^ %.4e km p^\n", r_inertial(1), r_inertial(2));
fprintf("v_inertial before maneuver: %.4e km/s e^ %.4e km/s p^\n", v_inertial_old(1), v_inertial_old(2));

% calc EA
EA_old = 2*atan(sqrt((1-e)/(1+e))*tand(TA_old/2));
fprintf("Eccentric Anomaly: %.4e rad\n", EA_old);
fprintf("Eccentric Anomaly: %.4e deg\n", rad2deg(EA_old));

%calc time since periapsis
t_since_tp = sqrt((a^3)/Gm_earth)*(EA_old - e*sin(EA_old));
fprintf("time since periapsis: %.4e sec\n", t_since_tp);
fprintf("time since periapsis: %.4e hours\n", t_since_tp/3600);
fprintf("time since periapsis: %.4e days\n", t_since_tp/3600/24);



% after maneuver properties
rp_new = 2*R_earth;
e_new = 0.4;

% calc new semi major axis
a_new = rp_new/(1-e_new);
fprintf("semi major axis after maneuver: %.4e\n", a_new);

% calculate semilatus rectum
p_new = a_new*(1-(e_new^2));
fprintf("semilatus rectum after maneuver: %.4e km\n", p_new);

%calc specific energy
specific_energy_new = -Gm_earth/(2*a_new);
fprintf("specific energy after maneuver: %.4e km^2/sec^2\n", specific_energy_new);

% calc v
v_new = sqrt(2*(specific_energy_new+(Gm_earth/r)));
fprintf("Velocity mag after maneuver: %.4e km/sec\n", v_new);

% calc new gamma
FPA_new = -acosd(sqrt(Gm_earth*p_new)/(r*v_new)); % choose negative
fprintf("FPA after maneuver: %.4e deg\n", FPA_new);

%calc d_FPA
d_FPA = FPA_new - FPA_old;
fprintf("delta FPA: %.4e deg\n", d_FPA);

% calc v_vector_new
v_vector_new = [ v_new*sind(FPA_new) v_new*cosd(FPA_new) ];
fprintf("v_vector after maneuver: %.4e km/sec r^ %.4e km.sec theta^\n", v_vector_new(1), v_vector_new(2));

% find v_new
d_v = sqrt((v_old^2) + (v_new^2) - (2*v_old*v_new*cosd(d_FPA)));
fprintf("delta v magnitude: %.4e km/sec\n", d_v);

%calc alpha
beta_angle = asind(v_new*sind(d_FPA)/d_v);
alpha_angle = 180-beta_angle;
fprintf("beta angle: %.4e deg\n", beta_angle);
fprintf("alpha angle: %.4e deg\n", alpha_angle);

% check beta and alpha using cosine rule
beta_check = acosd((v_new^2 - v_old^2 - d_v^2)/(-2*v_old*d_v));
alpha_check = 180-beta_check;
fprintf("beta check: %.4e deg\n", beta_check);
fprintf("alpha check: %.4e deg\n", alpha_check);

d_v_vector_vnb = [ d_v*cosd(alpha_check) d_v*sind(alpha_check)];
fprintf("d_v_vector in VNB: %.4e km/sec V^ %.4e km.sec B^\n", d_v_vector_vnb(1), d_v_vector_vnb(2));

%convert from vnb to polar:
TM_vnb_to_polar = [ sind(FPA_old) cosd(FPA_old); cosd(FPA_old) -sind(FPA_old)];
d_v_vector_polar = d_v_vector_vnb*TM_vnb_to_polar;
fprintf("d_v_vector in polar: %.4e km/sec r^ %.4e km.sec theta^\n", d_v_vector_polar(1), d_v_vector_polar(2));

% calc new eccentricity
e_check = sqrt(1-(p_new/a_new));
fprintf("eccentricity check after maneuver: %.4e\n", e_check);

% calc rp
rp_check = a_new*(1-e_check);
fprintf("rp check after maneuver: %.4e\n", rp_check);
rp_check_earth_ratio = rp_check/R_earth;
fprintf("rp check earth ratio after maneuver: %.4e\n", rp_check_earth_ratio);

% calc ra_new
ra_new = a_new*(1+e_new);
fprintf("ra after maneuver: %.4e\n", ra_new);

% calc new period
period_new = 2*pi*sqrt((a_new^3)/Gm_earth);
fprintf("period after maneuver: %.4e sec\n", period_new);
fprintf("period after maneuver: %.4e hours\n", period_new/3600);
fprintf("period after maneuver: %.4e days\n", period_new/3600/24);

% calc old TA
TA_new = -acosd((p_new/(r*e_new))-(1/e_new)); % negative because descending
fprintf("True Anomaly after maneuver: %.4e deg\n", TA_new);
fprintf("True Anomaly after maneuver: %.4e rad\n", deg2rad(TA_new));

% calc EA
EA_new = 2*atan(sqrt((1-e_new)/(1+e_new))*tand(TA_new/2));
fprintf("Eccentric Anomaly after maneuver: %.4e rad\n", EA_new);
fprintf("Eccentric Anomaly after maneuver: %.4e deg\n", rad2deg(EA_new));

%calc time since periapsis
t_since_tp_new = sqrt((a_new^3)/Gm_earth)*(EA_new - e_new*sin(EA_new));
fprintf("time since periapsis: %.4e sec\n", t_since_tp_new);
fprintf("time since periapsis: %.4e hours\n", t_since_tp_new/3600);
fprintf("time since periapsis: %.4e days\n", t_since_tp_new/3600/24);

d_little_omega = -TA_new + TA_old;
fprintf("shift in little omega: %.4e deg\n", d_little_omega);

% plot the original orbit
b = a * sqrt( 1 - (e^2) );
fprintf('old Semi-major Axis (a): %.7e[km]\n', a);
fprintf('old Semi-minor Axis (b): %.7e[km]\n', b);

% original orbit
t = linspace(0, 2*pi);
x = a * cos(t) - (a*e);
y = b * sin(t);
plot(x,y);
hold on;

% % new orbit
% b_new = a_new * sqrt( 1 - (e_new^2) );
% theta = deg2rad(d_little_omega);
% fprintf('new Semi-minor Axis (b_new): %.7e[km]\n', b_new);
% % use the general parametric equation for rotated ellipse
% x_new = a_new*cos(t)*cos(theta) - b_new*sin(t)*sin(theta) - (a_new*e_new*cos(theta));
% y_new = b_new*sin(t)*cos(theta) + a_new*cos(t)*sin(theta) - (a_new*e_new*sin(theta));
% plot(x_new,y_new, '-.');
% hold on;

% center lines
xline(0);
yline(0);
xline(-a*e)
hold on;

daspect([1 1 1]);
xlim([-3.5e4 1e4]);
ylim([-2e4 2e4]);
grid on;
hold off;

%legend("Before maneuver", "After maneuver")
title("Orbit of spacecraft around Earth-Lillian Shido")
xlabel("Distance [km]")
ylabel("Distance [km]")


