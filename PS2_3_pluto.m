% Given constants
Gm_pluto = 981.601;
Gm_charon = 119.480;
D_pluto = 1162;
D_charon = 606;
r_pluto_charon = 19596;
v_pluto = -0.025717;
v_charon = +0.211319;
G = 6.6743e-20;

% Calculated constants
m_pluto = Gm_pluto / G;
m_charon = Gm_charon / G;


% Problem 3a: CM



% Problem 3b: angular momentum
v_pluto_charon = v_charon - v_pluto;
fprintf("Pluto-Charon Velocity: %.5e\n", v_pluto_charon);

AV = v_pluto_charon / r_pluto_charon;
fprintf("Angular Velocity: %.5e\n", AV);


% Problem 3c: system linear momentum
LM_pluto = Gm_pluto*v_pluto/G;
LM_charon = Gm_charon*v_charon/G;

LM = LM_pluto + LM_charon;
fprintf("System Linear Momentum: %.5e\n", LM);

v_CM = LM * G / (Gm_pluto + Gm_charon);
fprintf("Velocity of CM: %.5e\n", v_CM);


% Problem 3d: Angular momentum constant C3
C3 = r_pluto_charon^2 * AV * (m_pluto*m_charon) / (m_pluto + m_charon);
fprintf("C3: %.5e\n", C3);


% Problem 3e: Energy constant C4
T = 1/2 * (v_pluto_charon^2) * (m_pluto * m_charon) / (m_pluto + m_charon);
fprintf("Kinetic Energy: %.5e\n", T);
U = G * m_pluto * m_charon / r_pluto_charon;
fprintf("Gravitational Potential: %.5e\n", U);
C4 = T - U;
fprintf("C4: %.5e\n", C4);
