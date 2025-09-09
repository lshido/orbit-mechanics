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

%solve f and g
f = 1 - (r_new/p*(1 - cosd(d_TA)));
g = r_new*r_1/sqrt(Gm_earth*p)*sind(d_TA);
f_dot = ((dot(r_inertial, v_inertial)/(p*r_1))*(1-cosd(d_TA))) - ((1/r_1)*sqrt(Gm_earth/p)*sind(d_TA));
g_dot = 1-((r_1/p)*(1-cosd(d_TA)));
fprintf("f: %.4e rad\n", f);
fprintf("g: %.4e sec\n", g);
fprintf("f_dot: %.4e 1/sec\n", f_dot);
fprintf("g_dot: %.4e rad \n", g_dot);

% solve f and g using EA
f_EA = 1-((a/r_1)*(1-cos(d_EA)));
g_EA = (d_t) - sqrt((a^3)/Gm_earth)*(d_EA - sin(d_EA));
f_dot_EA = -sqrt(Gm_earth*a)/(r_new*r_1)*sin(d_EA);
g_dot_EA = 1-(a/r_new)*(1-cos(d_EA));
fprintf("f_EA: %.4e rad\n", f_EA);
fprintf("g_EA: %.4e sec\n", g_EA);
fprintf("f_dot_EA: %.4e 1/sec\n", f_dot_EA);
fprintf("g_dot_EA: %.4e rad \n", g_dot_EA);

r_inertial_2 = f*r_inertial + g*v_inertial;
v_inertial_2 = f_dot*r_inertial + g_dot*v_inertial;
fprintf("r_inertial at t2: %.4e km x^ %.4e km y^ %.4e km z^\n", r_inertial_2(1), r_inertial_2(2), r_inertial_2(3));
fprintf("v_inertial at t2: %.4e km/s x^ %.4e km/s y^ %.4e km/s z^\n", v_inertial_2(1), v_inertial_2(2), v_inertial_2(3));

r_polar_1 = [ r_1 0 0];
v_polar_1 = [-r_dot r_theta_dot 0 ];
f_dot_polar = ((dot(r_polar_1, v_polar_1)/(p*r_1))*(1-cosd(d_TA))) - ((1/r_1)*sqrt(Gm_earth/p)*sind(d_TA));
g_dot_polar = 1-((r_1/p)*(1-cosd(d_TA)));
v_polar_2 = f_dot_polar*r_polar_1 + g_dot_polar*v_polar_1;
fprintf("v_polar at t2: %.4e km/s r^ %.4e km/s theta^ %.4e km/s h^\n", v_polar_2(1), v_polar_2(2), v_polar_2(3));

% calc r_dot 
r_dot_2 = sqrt((v_new^2)-(p*Gm_earth/(r_new^2)));
fprintf("r_dot_2 mag at t1: %.4e km/sec\n", r_dot_2);

% calc r_theta_dot
r_theta_dot_2 = sqrt(p*Gm_earth)/r_new;
fprintf("r_theta_dot_2 mag at t1: %.4e km/sec\n", r_theta_dot_2);

% calc FPA
FPA_2 = atand(r_dot_2/r_theta_dot_2);
fprintf("FPA at t2: %.4e deg\n", FPA_2);
fprintf("FPA_GMAT at t2: %.4e deg\n", 90 - FPA_2);

% Plot the orbit
b = a * sqrt( 1 - (e^2) );
fprintf('Semi-major Axis (a): %.7e[km]\n', a);
fprintf('Semi-minor Axis (b): %.7e[km]\n', b);

r_p = a * (1-e);
fprintf('radius at periapsis: %.7e[km]\n', r_p);

% ellipse
t = linspace(0, 2*pi);
x = a * cos(t) - (a*e);
y = b * sin(t);
plot(x,y);
hold on;

% reference cirle
x_c = a * cos(t) - (a*e);
y_c = a * sin(t);
plot(x_c,y_c,'--');

% center lines
xline(0);
yline(0);


axis equal;
grid on;
hold off;
axis auto;
title("Orbit of spacecraft around Earth-Lillian Shido")
xlabel("Distance [km]")
ylabel("Distance [km]")