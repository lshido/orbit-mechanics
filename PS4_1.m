R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
e = 0.6;
p = 6 * R_earth;
theta_star_0 = 90; % [deg]

% 1(a)
fprintf('Problem 1a:\n')
% Find semi-major axis a:
a = p / (1- (e^2) );
fprintf('Semi-major axis (a): %.7e[km]\n', a);
fprintf('a with respect to R_earth: %.7e\n', a/R_earth);

% Find r at periapsis r_p:
r_p = a * (1-e);
fprintf('Radius at periapsis (r_p)): %.7e[km]\n', r_p);
fprintf('With respect to R_earth: %.7e\n', r_p/R_earth);

% Find r at apoapsis r_a:
r_a = a * (1+e);
fprintf('Radius at apoapsis (r_a)): %.7e[km]\n', r_a);
fprintf('With respect to R_earth: %.7e\n', r_a/R_earth);

% Find period P:
P = 2 * pi * sqrt( (a^3)/(Gm_earth) );
fprintf('Period (P))): %.7e[sec]\n', P);
fprintf('Period (P))): %.7e[hr]\n', P/3600);

% Find specific energy curly_E:
curly_E = - Gm_earth / (2 * a);
fprintf('Specific Energy (curly_E): %.7e[km^2/s^2]\n', curly_E);

% Find r_0:
r_0 = p / ( 1 + ( e * cosd(theta_star_0)));
fprintf('Initial r (r_0): %.7e[km]\n', r_0);

% Find v_0:
v_0 = sqrt( 2 * (curly_E + (Gm_earth/r_0)));
fprintf('Initial v (v_0): %.7e[km/sec]\n', v_0);

% Find eccentric_anomaly:
EA_0 = acosd( (a-r_0) / (a*e) );
fprintf('Initial Eccentric Anomaly (EA_0): %.7e[deg]\n', EA_0);

% Find scalar r_dot:
scalar_r_dot = sqrt( (v_0^2) - ( (p * Gm_earth) / (r_0^2) ) );
fprintf('Scalar r dot (r_dot): %.7e[km/sec]\n', scalar_r_dot);

% Find scalar r_times_theta_dot:
r_times_theta_dot = sqrt(p * Gm_earth)/r_0;
fprintf('Scalar r times theta dot (r_times_theta_dot): %.7e[km/sec]\n', r_times_theta_dot);

% Find flight path angle:
FPA = atand( scalar_r_dot / r_times_theta_dot );
fprintf('Flight Path Angle (FPA): %.7e[deg]\n', FPA);

% Find circular velocity:
v_c = sqrt( Gm_earth / r_0 );
fprintf('Circular Velocity (v_c): %.7e[km/sec]\n', v_c);

% Find escape speed:
v_escape = sqrt(2) * v_c;
fprintf('Escape Speed (v_escape): %.7e[km/sec]\n', v_escape);

% 1(b)
fprintf('\n')
fprintf('Problem 1b:\n')
EA_f = 225; % [deg]

% Find final r:
r_f = a * ( 1 - ( e * cosd(EA_f) ) );
fprintf('Final position (r_f): %.7e[km]\n', r_f);

% Find v_f:
v_f = sqrt( 2 * (curly_E + (Gm_earth/r_f)));
fprintf('Final velocity (v_f): %.7e[km/sec]\n', v_f);

% Find final true anomaly:
TA_f = -acosd( (p / (e * r_f)) - (1 / e) );
fprintf('Final True Anomaly (TA_f): %.7e[deg]\n', TA_f);

% Find scalar r_dot_f:
scalar_r_dot_f = sqrt( (v_f^2) - ( (p * Gm_earth) / (r_f^2) ) );
fprintf('Final Scalar r dot (r_dot_f): %.7e[km/sec]\n', scalar_r_dot_f);

% Find final scalar r_times_theta_dot:
r_times_theta_dot_f = sqrt(p * Gm_earth)/r_f;
fprintf('Final Scalar r times theta dot (r_times_theta_dot_f): %.7e[km/sec]\n', r_times_theta_dot_f);

% Find final flight path angle:
FPA_f = atand( -scalar_r_dot_f / r_times_theta_dot_f );
fprintf('Final Flight Path Angle (FPA_f): %.7e[deg]\n', FPA_f);

% Find position vector w.r.t. e_hat and p_hat:
r_f_e_p = r_f * [cosd(TA_f) sind(TA_f)];
fprintf('Position in e_hat and p_hat @ t_f: %.7e[km]e_hat %.7e[km]p_hat\n', r_f_e_p(1), r_f_e_p(2));

% 1(c)
fprintf('\n')
fprintf('Problem 1c:\n')
% Find t_0 - t_p:
t_0_from_t_p = ( EA_0 - (e * sind(EA_0))) / sqrt(Gm_earth / (a^3));
fprintf('t_0 relative to t_p (sec): %.7e[sec]\n', t_0_from_t_p);
fprintf('t_0 relative to t_p (hours): %.7e[hours]\n', t_0_from_t_p/3600);
fprintf('t_0 relative to t_p (days): %.7e[days]\n', t_0_from_t_p/3600/24);

% Find t_f - t_p:
t_f_from_t_p = ( EA_f - (e * sind(EA_f))) / sqrt(Gm_earth / (a^3));
fprintf('t_f relative to t_p (sec): %.7e[sec]\n', t_f_from_t_p);
fprintf('t_f relative to t_p (hours): %.7e[hours]\n', t_f_from_t_p/3600);
fprintf('t_f relative to t_p (days): %.7e[days]\n', t_f_from_t_p/3600/24);

% Find t_f - t_0:
t_f_from_t_0 = t_f_from_t_p - t_0_from_t_p;
fprintf('t_f - t_0 (sec): %.7e[sec]\n', t_f_from_t_0);
fprintf('t_f - t_0 (hours): %.7e[hours]\n', t_f_from_t_0/3600);
fprintf('t_f - t_0 (days): %.7e[days]\n', t_f_from_t_0/3600/24);

% Find delta True anomaly:
d_TA = (360+TA_f) - theta_star_0;
fprintf('Delta True Anomaly (d_TA): %.7e[deg]\n', d_TA);

% Find delta Eccentric anomaly:
d_EA = EA_f - EA_0;
fprintf('Delta Eccentric Anomaly (d_EA): %.7e[deg]\n', d_EA);