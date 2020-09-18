% Givens
m_sc = 600;
Gm_earth = 398600.4415;
G = 6.6743e-20;
r_earth = 6378.1363;

% Givens at t=0
altitude = 8560;
v_I = [2.11 4.89 0];

% Calculated constants
m_earth = Gm_earth / G;

% Problem 3a
r = [r_earth + altitude 0 0];
h = cross(r, v_I);
fprintf("Specific Angular Momentum (h): %.2er %.2etheta %.2eh\n",...
    h(1), h(2), h(3));

C_scaling = (m_earth + m_sc) / (m_earth * m_sc);
fprintf("Scaling Constant: %.2e\n", C_scaling);

C3 = h / C_scaling;
fprintf("Angular Momentum Constant (C3): %.2er %.2etheta %.2eh\n",...
    C3(1), C3(2), C3(3));

T = 1/2 / C_scaling * dot(v_I, v_I);
U = - Gm_earth * m_sc / norm(r);
C4 = T - U;
fprintf("T: %.2e\n", T);
fprintf("U: %.2e\n", U);
fprintf("Energy Constant (C4): %.2e\n", C4);

E = C_scaling * C4;
fprintf("Specific Energy (E): %.2e\n", E);

A_dot = norm(h) / 2;
fprintf("Areal Velocity (A_dot): %.2e\n", A_dot);

% Problem 3c
U_dash = 1/2 * norm(v_I)^2 - E;
mu = U_dash * norm(r);
fprintf("mu: %.2e\n", mu);

e_vector = (cross(v_I, h) / mu) - (r / norm(r));
fprintf("Eccentricity (e): %.2e\n", norm(e_vector));

p = (norm(h))^2 / mu;
fprintf("Semilatus Rectum (p): %.2e\n", p);

a = p / (1 - (norm(e))^2);
fprintf("Semimajor Axis (a): %.2e\n", a);

period = 2 * pi * sqrt(a^3 / mu);
fprintf("Period (P): %.2e\n", period);

theta_star = acosd((p/norm(r))/norm(e));
fprintf("True Anomaly (theta_star): %.2e\n", theta_star);
