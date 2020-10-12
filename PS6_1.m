R_mars = 3397; % km
Gm_mars = 42828.314258067; % km^3/sec^2

r_p = 1.1 * R_mars;
r_a = 6.0 * R_mars;
TA_c = 90; % deg


% calc e and a
syms e a
equations = [ a*(1-e) == r_p, a*(1+e) == r_a ];
variables = [ a e ];
[solved_a, solved_e] = solve(equations, variables);
fprintf("semi-major axis: %.4e km\n", solved_a);
fprintf("eccentricity: %.4e\n", solved_e);

% calculate semilatus rectum
p = solved_a*(1-(solved_e^2));
fprintf("semilatus rectum: %.4e km\n", p);

% calc TA at r
r_1 = 4.5 * R_mars;
fprintf("Radius at t1: %.4e km\n", r_1);
TA_1 = acosd((p/(r_1*solved_e))-(1/solved_e));
fprintf("True Anomaly at t1: %.4e deg\n", TA_1);
fprintf("True Anomaly at t1: %.4e rad\n", deg2rad(TA_1));

%calc energy
specific_energy = -Gm_mars/(2*solved_a);
fprintf("specific energy: %.4e km^2/sec^2\n", specific_energy);

% calc v
v_1 = sqrt(2*(specific_energy+(Gm_mars/r_1)));
fprintf("Velocity mag at M: %.4e km/sec\n", v_1);

% calc r_dot 
r_dot = sqrt((v_1^2)-(p*Gm_mars/(r_1^2)));
fprintf("r_dot mag at M: %.4e km/sec\n", r_dot);

% calc r_theta_dot
r_theta_dot = sqrt(p*Gm_mars)/r_1;
fprintf("r_theta_dot mag at M: %.4e km/sec\n", r_theta_dot);

% calc FPA
FPA = atand(r_dot/r_theta_dot);
fprintf("FPA at M: %.4e deg\n", FPA);

% calc EA
EA_1 = 2*atan(sqrt((1-solved_e)/(1+solved_e))*tand(TA_1/2));
fprintf("Eccentric Anomaly: %.4e rad\n", EA_1);
fprintf("Eccentric Anomaly: %.4e deg\n", rad2deg(EA_1));

% calc MA
MA_1 = EA_1 - solved_e*sin(EA_1);
fprintf("Mean Anomaly at t1: %.4e rad\n", MA_1);
fprintf("Mean Anomaly at t1: %.4e deg\n", rad2deg(MA_1));

%calc time since periapsis
t_since_t_p = MA_1*sqrt((solved_a^3)/Gm_mars);
fprintf("time since periapsis: %.4e sec\n", t_since_t_p);
fprintf("time since periapsis: %.4e hours\n", t_since_t_p/3600);
fprintf("time since periapsis: %.4e days\n", t_since_t_p/3600/24);

% calc EA at current
EA_c = 2*atan(sqrt((1-solved_e)/(1+solved_e))*tand(TA_c/2));
fprintf("Eccentric Anomaly at current: %.4e rad\n", EA_c);
fprintf("Eccentric Anomaly at current: %.4e deg\n", rad2deg(EA_c));

% calc MA at current
MA_c = EA_c - solved_e*sin(EA_c);
fprintf("Mean Anomaly at current: %.4e rad\n", MA_c);
fprintf("Mean Anomaly at current: %.4e deg\n", rad2deg(MA_c));

%calc time since periapsis
t_c_since_t_p = MA_c*sqrt((solved_a^3)/Gm_mars);
fprintf("time since periapsis: %.4e sec\n", t_c_since_t_p);
fprintf("time since periapsis: %.4e hours\n", t_c_since_t_p/3600);
fprintf("time since periapsis: %.4e days\n", t_c_since_t_p/3600/24);

wait_time = t_since_t_p - t_c_since_t_p;
fprintf("wait time: %.4e sec\n", wait_time);
fprintf("wait time: %.4e hours\n", wait_time/3600);
fprintf("wait time: %.4e days\n", wait_time/3600/24);

% calc new v
r_new = r_1;
a_new = r_1;
v_new = sqrt(2*((Gm_mars/r_new) - (Gm_mars/(2*a_new))));
fprintf("velocity magnitude after maneuver: %.4e km/sec\n", v_new);

% calc new gamma
p_new = r_new;
gamma_new = acosd(sqrt(p_new*Gm_mars)/(r_new*v_new));
fprintf("FPA after maneuver: %.4e deg\n", gamma_new);

% calc delta v
d_gamma = FPA;
d_v = sqrt((v_new^2) + (v_1^2) - (2*v_new*v_1*cosd(d_gamma)));
fprintf("delta v magnitude: %.4e km/sec\n", d_v);

%calc alpha
beta_angle = asind(v_new*sind(d_gamma)/d_v);
fprintf("beta angle: %.4e deg\n", beta_angle);
fprintf("alpha angle: %.4e deg\n", 180-beta_angle);

% Plot the orbit
b = solved_a * sqrt( 1 - (solved_e^2) );
fprintf('Semi-major Axis (a): %.7e[km]\n', solved_a);
fprintf('Semi-minor Axis (b): %.7e[km]\n', b);

% ellipse
t = linspace(0, 2*pi);
x = solved_a * cos(t) - (solved_a*solved_e);
y = b * sin(t);
plot(x,y);
hold on;

% new orbit
x_new = r_new * cos(t);
y_new = r_new * sin(t);
plot(x_new,y_new, '-.');
hold on;

% center lines
xline(0);
yline(0);
hold on;

axis equal;
grid on;
hold off;
axis auto;
title("Orbit of spacecraft around Mars-Lillian Shido")
xlabel("Distance [km]")
ylabel("Distance [km]")

