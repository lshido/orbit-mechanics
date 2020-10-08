R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

a = 20*R_earth;
fprintf('Semi-major axis: %.4e[km]\n', a);
e = 0.6;
inclination = 34; %[deg]
big_omega = 45; %[deg]
small_omega = 30; %[deg]
theta = 235; %[deg]

% calc TA
TA = theta - small_omega;
fprintf('True Anomaly: %.4e[deg]\n', TA);
fprintf('True Anomaly: %.4e[deg]\n', TA-360);

% calc semilatus rectum 
p = a * (1 - e^2);
fprintf('Semi-latus rectum: %.4e[km]\n', p);

r_1 = p/(1+e*cosd(TA));
fprintf('r_1: %.4e[km]\n', r_1);

%calc energy
specific_energy = -Gm_earth/(2*a);
fprintf("specific energy: %.4e km^2/sec^2\n", specific_energy);

% calc v_1
v_1 = sqrt(2*(specific_energy+(Gm_earth/r_1)));
fprintf("Velocity mag at t1: %.4e km/sec\n", v_1);

%calc period
period = 2*pi*sqrt((a^3)/Gm_earth);
fprintf("period: %.4e sec\n", period);
fprintf("period: %.4e hours\n", period/3600);
fprintf("period: %.4e days\n", period/3600/24);

% calc r_dot 
r_dot = sqrt((v_1^2)-(p*Gm_earth/(r_1^2)));
fprintf("r_dot mag at t1: %.4e km/sec\n", r_dot);

% calc r_theta_dot
r_theta_dot = sqrt(p*Gm_earth)/r_1;
fprintf("r_theta_dot mag at t1: %.4e km/sec\n", r_theta_dot);

% calc FPA
FPA = atand(-r_dot/r_theta_dot);
fprintf("FPA at t1: %.4e deg\n", FPA);

% calc EA
EA = 2*atan(sqrt((1-e)/(1+e))*tand(TA/2));
fprintf("Eccentric Anomaly: %.4e rad\n", EA);
fprintf("Eccentric Anomaly: %.4e deg\n", rad2deg(EA));

% calc MA
MA = EA - e*sin(EA);
fprintf("Mean Anomaly at t1: %.4e rad\n", MA);
fprintf("Mean Anomaly at t1: %.4e deg\n", rad2deg(MA));

%calc time since periapsis
t_since_t_p = MA*sqrt((a^3)/Gm_earth);
fprintf("time since periapsis: %.4e sec\n", t_since_t_p);
fprintf("time since periapsis: %.4e hours\n", t_since_t_p/3600);
fprintf("time since periapsis: %.4e days\n", t_since_t_p/3600/24);

% express r and v in n_hat
TM_1 = [ cosd(theta) sind(theta) 0; -sind(theta) cosd(theta) 0; 0 0 1 ];
TM_2 = [ 1 0 0; 0 cosd(inclination) sind(inclination); 0 -sind(inclination) cosd(inclination) ];
r_qhat = [r_1 0 0] * TM_1;
r_nhat = r_qhat * TM_2;
fprintf("r_nhat at t1: %.4e km nx^ %.4e km ny^ %.4e km nz^\n", r_nhat(1), r_nhat(2), r_nhat(3));

v_qhat = [-r_dot r_theta_dot 0] * TM_1;
v_nhat = v_qhat * TM_2;
fprintf("v_nhat at t1: %.4e km/s nx^ %.4e km/s ny^ %.4e km/s nz^\n", v_nhat(1), v_nhat(2), v_nhat(3));

% express r and v in inertial units
TM_3 = [ cosd(big_omega) sind(big_omega) 0; -sind(big_omega) cosd(big_omega) 0; 0 0 1 ];
r_inertial = r_nhat * TM_3;
v_inertial = v_nhat * TM_3;
fprintf("r_inertial at t1: %.4e km x^ %.4e km y^ %.4e km z^\n", r_inertial(1), r_inertial(2), r_inertial(3));
fprintf("v_inertial at t1: %.4e km/s x^ %.4e km/s y^ %.4e km/s z^\n", v_inertial(1), v_inertial(2), v_inertial(3));

% calc 3 days later
d_t = 3*24*3600;
new_t_since_t_p = t_since_t_p + d_t;
fprintf("new time since periapsis: %.4e sec\n", new_t_since_t_p);
fprintf("new time since periapsis: %.4e hours\n", new_t_since_t_p/3600);
fprintf("new time since periapsis: %.4e days\n", new_t_since_t_p/3600/24);

%calc M
M_new = sqrt(Gm_earth/(a^3))*new_t_since_t_p;
fprintf("M_new: %.4e rads\n", M_new);
fprintf("M_new: %.4e deg\n", rad2deg(M_new));

%solve for EA
syms EA_new
Kepler_equation_new = EA_new - e*sin(EA_new) == M_new;
solved_EA_new = vpasolve(Kepler_equation_new, EA_new);
fprintf("Eccentric Anomaly: %.4e rads\n", solved_EA_new);
fprintf("Eccentric Anomaly: %.4e degrees\n", rad2deg(solved_EA_new));

%calc TA
%solve for new TA
TA_new = 2*atan(sqrt((1+e)/(1-e))*tan(solved_EA_new/2));
fprintf("True Anomaly: %.4e rad\n", TA_new);
fprintf("True Anomaly: %.4e deg\n", rad2deg(TA_new));

%calc new r
r_new = p/(1+(e*cos(TA_new)));
fprintf("r_new: %.4e km\n", r_new);

% calc new v
v_new = sqrt(2*(specific_energy+(Gm_earth/r_new)));
fprintf("Velocity mag at t2: %.4e km/sec\n", v_new);

d_TA = rad2deg(TA_new) - TA;
fprintf("Difference in True Anomaly: %.4e deg\n", d_TA);
fprintf("Difference in True Anomaly: %.4e rad\n", deg2rad(d_TA));

d_EA = solved_EA_new - EA;
fprintf("Difference in Eccentric Anomaly: %.4e rad\n", d_EA);
fprintf("Difference in Eccentric Anomaly: %.4e deg\n", rad2deg(d_EA));
