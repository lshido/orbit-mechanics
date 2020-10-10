R_mars = 3397; % km
Gm_mars = 42828.314258067; % km^3/sec^2

r_p = 1.5 * R_mars;
r_a = 6.5 * R_mars;
M = -90; %deg
fprintf("radius at periapsis: %.4e km\n", r_p);
fprintf("radius at apoapsis: %.4e km\n", r_a);

syms e a
equations = [a*(1-e) == r_p, a*(1+e) == r_a];
variables = [a e];
[solved_a, solved_e] = solve(equations, variables);
fprintf("semi-major axis: %.4e km\n", solved_a);
fprintf("eccentricity: %.4e\n", solved_e);

% calculate semilatus rectum
p = solved_a*(1-(solved_e^2));
fprintf("semilatus rectum: %.4e km\n", p);

% calculate specific angular momentum
h = sqrt(p * Gm_mars);
fprintf("specific angular momentum: %.4e km^2/sec\n", h);

%calc period
period = 2*pi*sqrt((solved_a^3)/Gm_mars);
fprintf("period: %.4e sec\n", period);
fprintf("period: %.4e hours\n", period/3600);

%calc energy
specific_energy = -Gm_mars/(2*solved_a);
fprintf("specific energy: %.4e km^2/sec^2\n", specific_energy);

%calc time since periapsis
t_since_t_p = deg2rad(M)*sqrt((solved_a^3)/Gm_mars);
fprintf("time since periapsis: %.4e sec\n", t_since_t_p);
fprintf("time since periapsis: %.4e hours\n", t_since_t_p/3600);

%solve for EA - Lillian Shido
M = -90; %deg
syms EA
Kepler_equation = EA - solved_e*sin(EA) == deg2rad(M);
solved_EA = vpasolve(Kepler_equation, EA);
fprintf("Eccentric Anomaly: %.4e rads\n", solved_EA);
fprintf("Eccentric Anomaly: %.4e degrees\n", rad2deg(solved_EA));

% calc TA
TA = 2*atan(sqrt((1+solved_e)/(1-solved_e))*tan(solved_EA/2));
fprintf("True Anomaly: %.4e rad\n", TA);
fprintf("True Anomaly: %.4e deg\n", rad2deg(TA));

%calc r
r = p/(1+(solved_e*cos(TA)));
fprintf("Radius mag at M: %.4e km\n", r);

% calc v
v = sqrt(2*(specific_energy+(Gm_mars/r)));
fprintf("Velocity mag at M: %.4e km/sec\n", v);

% calc r_dot 
r_dot = sqrt((v^2)-(p*Gm_mars/(r^2)));
fprintf("r_dot mag at M: %.4e km/sec\n", r_dot);

% calc r_theta_dot
r_theta_dot = sqrt(p*Gm_mars)/r;
fprintf("r_theta_dot mag at M: %.4e km/sec\n", r_theta_dot);

% calc FPA
FPA = atand(-r_dot/r_theta_dot);
fprintf("FPA at M: %.4e deg\n", FPA);

TM = [ cos(TA) sin(TA); -sin(TA) cos(TA) ];
r_inertial = [r 0] * TM;
fprintf("r_inertial at M: %.4e km e^ %.4e km p^\n", r_inertial(1), r_inertial(2));

v_inertial = [-r_dot r_theta_dot] * TM;
fprintf("v_inertial at M: %.4e km/s e^ %.4e km/s p^\n", v_inertial(1), v_inertial(2));

new_t_since_t_p = t_since_t_p + (0.5 * period);
fprintf("new time since periapsis: %.4e sec\n", new_t_since_t_p);
fprintf("new time since periapsis: %.4e hours\n", new_t_since_t_p/3600);

%solve for EA - Lillian Shido
M_new = sqrt(Gm_mars/(solved_a^3))*new_t_since_t_p;
fprintf("M_new: %.4e rads\n", M_new);
fprintf("M_new: %.4e deg\n", rad2deg(M_new));
syms EA_new
Kepler_equation_new = EA_new - solved_e*sin(EA_new) == M_new;
solved_EA_new = vpasolve(Kepler_equation_new, EA_new);
fprintf("Eccentric Anomaly: %.4e rads\n", solved_EA_new);
fprintf("Eccentric Anomaly: %.4e degrees\n", rad2deg(solved_EA_new));

%solve for new TA
TA_new = 2*atan(sqrt((1+solved_e)/(1-solved_e))*tan(solved_EA_new/2));
fprintf("True Anomaly: %.4e rad\n", TA_new);
fprintf("True Anomaly: %.4e deg\n", rad2deg(TA_new));

%calc new r
r_new = p/(1+(solved_e*cos(TA_new)));
fprintf("r_new: %.4e km\n", r_new);

%solve f and g
d_TA = TA_new - TA;
fprintf("Difference in True Anomaly: %.4e deg\n", rad2deg(d_TA));
f = 1 - (r_new/p*(1 - cos(d_TA)));
g = r_new*r/sqrt(Gm_mars*p)*sin(d_TA);
f_dot = ((dot(r_inertial, v_inertial)/(p*r))*(1-cos(d_TA))) - ((1/r)*sqrt(Gm_mars/p)*sin(d_TA));
g_dot = 1-((r/p)*(1-cos(d_TA)));
fprintf("f: %.4e rad\n", f);
fprintf("g: %.4e sec\n", g);
fprintf("f_dot: %.4e 1/sec\n", f_dot);
fprintf("g_dot: %.4e rad \n", g_dot);

% solve f and g using EA
d_EA = solved_EA_new - solved_EA;
fprintf("Difference in Eccentric Anomaly: %.4e deg\n", rad2deg(d_EA));
f_EA = 1-((solved_a/r)*(1-cos(d_EA)));
g_EA = (0.5*period) - sqrt((solved_a^3)/Gm_mars)*(d_EA - sin(d_EA));
fprintf("f_EA: %.4e rad\n", f_EA);
fprintf("g_EA: %.4e sec\n", g_EA);

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

% reference cirle
x_c = solved_a * cos(t) - (solved_a*solved_e);
y_c = solved_a * sin(t);
plot(x_c,y_c,'--');

% center lines
xline(0);
yline(0);


axis equal;
grid on;
hold off;
axis auto;
title("Orbit of spacecraft around Mars-Lillian Shido")
xlabel("Distance [km]")
ylabel("Distance [km]")

% iterative solver for EA - Lillian Shido
M = 90; % [deg]
e = 0.625;
MA_error = 1;
EA_guess = deg2rad(M);
stepsize = 0.001;

while MA_error > 0.001
    if MA_error > 0
        EA_guess = EA_guess - stepsize;
    else
        EA_guess = EA_guess + stepsize;
    end
    MA_calced = EA_guess - solved_e*sin(EA_guess);
    MA_error = MA_calced - deg2rad(M);
    fprintf("Guessed EA: %.4e rad    MA_error: %.4e\n", EA_guess, MA_error);
end
fprintf("Calculated EA: %.4e rad\n", EA_guess);
fprintf("Calculated EA: %.4e degrees\n", rad2deg(EA_guess));
