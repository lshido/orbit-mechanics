R_moon = 1738.2;
Gm_moon = 4902.8005821478;


a = 7050;
r_p = R_moon + 800;
fprintf("r_p: %.4e km\n", r_p);

% specific energy
energy = Gm_moon/(2*a);
fprintf("specific energy: %.4e km^2/sec^2\n", energy);

% velocity at periapsis
v_p = sqrt(2*(energy + (Gm_moon/r_p)));
fprintf("v_p: %.4e km/sec\n", v_p);

% eccentricity
e = r_p/a + 1;
fprintf("eccentricity: %.4e \n", e);

% specific angular momentum
p = a*((e^2) - 1);
h = sqrt(p* Gm_moon);
fprintf("specific angular momentum: %.4e km^2/sec\n", h);

% flyby angle
delta = 2*asind(1/e);
fprintf("flyby angle: %.4e degrees\n", delta);

% b
theta_star_inf = acosd(-1/e);
b = (r_p + a)*sind(180 - theta_star_inf);
fprintf("b: %.4e km\n", b);

% v_infinity
v_inf = sqrt(2*energy);
fprintf("velocity at infinity: %.4e km/sec\n", v_inf);

% r at -60
theta_60 = -60;
r_60 = p/(1 + (e*cosd(theta_60)));
fprintf("r at -60: %.4e km\n", r_60);

% velocity at -60
v_60 = sqrt(2*(energy + (Gm_moon/r_60)));
fprintf("v at -60: %.4e km/sec\n", v_60);

% flight path angle gamma at -60
r_dot_60 = sqrt(v_60^2 - ((p*Gm_moon)/r_60^2));
r_theta_dot_60 = sqrt(p*Gm_moon) / r_60;
gamma_60 = atand(-r_dot_60/r_theta_dot_60);
fprintf("r_dot_60: %.4e km/sec\n", r_dot_60);
fprintf("r_theta_dot_60: %.4e km/sec\n", r_theta_dot_60);
fprintf("FPA: %.4e deg\n", gamma_60);

% Hyperbolic anomaly
H_60 = 2*atanh(sqrt((e-1)/(e+1))*tand(theta_60/2));
fprintf("H: %.4e deg\n", H_60);

% time till perilune

time_to_perilune = sqrt((a^3)/Gm_moon)*((e*sinh(H))-H);
fprintf("time to perilune: %.4e sec\n", time_to_perilune);
fprintf("time to perilune: %.4e hours\n", time_to_perilune/3600);

% hyperbola
t = linspace(deg2rad(-100), deg2rad(100));
x = a * -cosh(t) + (a+r_p);
y = p * sinh(t);
plot(x,y);
hold on;

% center lines
xline(a+r_p);
yline(0);

axis equal;
grid on;
hold off;
axis auto;
title("Hyperbolic Orbit of spacecraft around Earth")
xlabel("Distance [km]")
ylabel("Distance [km]")

% r at 100
theta_100 = 100;
r_100 = p/(1 + (e*cosd(theta_100)));
fprintf("r at 100: %.4e km\n", r_100);

% velocity at -60
v_100 = sqrt(2*(energy + (Gm_moon/r_100)));
fprintf("v at 100: %.4e km/sec\n", v_100);

% flight path angle gamma at -60
r_dot_100 = sqrt(v_100^2 - ((p*Gm_moon)/r_100^2));
r_theta_dot_100 = sqrt(p*Gm_moon) / r_100;
gamma_100 = atand(r_dot/r_theta_dot);
fprintf("r_dot_100: %.4e deg\n", r_dot_100);
fprintf("r_theta_dot_100: %.4e deg\n", r_theta_dot_100);
fprintf("FPA: %.4e deg\n", gamma_100);
