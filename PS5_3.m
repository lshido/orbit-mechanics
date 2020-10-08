R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;

a = 20*R_earth;
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
