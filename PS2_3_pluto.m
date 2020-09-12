% Given constants
Gm_pluto = 981.601;
Gm_charon = 119.480;
D_pluto = 1162;
D_charon = 606;
r_pluto_charon = [19596 0 0];
v_pluto = [0 -0.025717 0];
v_charon = [0 +0.211319 0];
G = 6.6743e-20;

% Calculated constants
m_pluto = Gm_pluto / G;
m_charon = Gm_charon / G;


% Problem 3a: CM
CM = (Gm_charon * r_pluto_charon) / (Gm_pluto + Gm_charon);
fprintf("Center of Mass: %.4ei %.4ej %.4ek\n", CM(1), CM(2), CM(3));


% Problem 3b: angular velocity
r_theta_dot = norm(r_pluto_charon) * norm(AV);
fprintf("r_theta_dot: %.4e\n", r_theta_dot);


v_pluto_charon = v_charon - v_pluto;
fprintf("Pluto-Charon Velocity: %.4ei %.4ej %.4ek\n", v_pluto_charon(1), v_pluto_charon(2), v_pluto_charon(3));

AV = v_pluto_charon / norm(r_pluto_charon);
fprintf("Angular Velocity: %.4ei %.4ej %.4ek\n", AV(1), AV(2), AV(3));


% Problem 3c: system linear momentum
LM_pluto = m_pluto*v_pluto;
fprintf("LM_pluto: %.4ei %.4ej %.4ek\n", LM_pluto(1), LM_pluto(2), LM_pluto(3));
LM_charon = m_charon*v_charon;
fprintf("LM_charon: %.4ei %.4ej %.4ek\n", LM_charon(1), LM_charon(2), LM_charon(3));

LM = LM_pluto + LM_charon;
fprintf("System Linear Momentum: %.4ei %.4ej %.4ek\n", LM(1), LM(2), LM(3));

v_CM = LM / (Gm_pluto + Gm_charon) * G;
fprintf("Velocity of CM: %.4e i %.4e j %.4e k\n", v_CM(1), v_CM(2), v_CM(3));


% Problem 3d: Angular momentum constant C3
C3_AV = norm(r_pluto_charon)^2 * AV * (m_pluto*m_charon) / (m_pluto + m_charon);
fprintf("C3 from angular velocity: %.4e i %.4e j %.4e k\n", C3_AV(1), C3_AV(2), C3_AV(3));

C3 = (m_pluto*m_charon) / (m_pluto + m_charon) * (cross(r_pluto_charon, v_pluto_charon));
fprintf("C3 from cross_product: %.4e i %.4e j %.4e k\n", C3(1), C3(2), C3(3));


% Problem 3e: Energy constant C4
T = 1/2 * dot(v_pluto_charon, v_pluto_charon) * (m_pluto * m_charon) / (m_pluto + m_charon);
fprintf("Kinetic Energy: %.4e\n", T);
U = G * m_pluto * m_charon / norm(r_pluto_charon);
fprintf("Gravitational Potential: %.4e\n", U);
C4 = T - U;
fprintf("C4: %.4e\n", C4);
