R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
r_earth_moon = 384400;
Gm_moon = 4902.8005821478;



D = [2, 10, 75, 200, 800, 1e12];

% Calculate semilatus rectum
p = 2 * (R_earth + 225);
fprintf('Semilatus Rectum: %.7e km\n', p);

% Calculate velocity
for i = 1:length(D);
    r = D(i) * R_earth;
    v = sqrt( 2*Gm_earth / r );
    r_dot = sqrt(v^2 - ((p*Gm_earth) / r^2) );
    r_theta_dot = sqrt(p*Gm_earth) / r;
    fprintf('Velocity at %d*R: %.7e km ( %.7e r^ %.7etheta^)\n', D(i), v, r_dot, r_theta_dot);
end;

% Calculate true anomaly
for i = 1:length(D);
    r = D(i) * R_earth;
    TA = acosd( ( 2*(R_earth + 225) / r) - 1 );
    fprintf('TA at %d*R: %.7e deg\n', D(i), TA);
    time_since_periapsis = (1/6)*sqrt(p^3/Gm_earth)*( (tand(TA/2)^3) + ( 3*tand(TA/2) ) );
    fprintf('t- t_p at %d*R: %.7e sec\n', D(i), time_since_periapsis);
    fprintf('t- t_p at %d*R: %.7e hours\n', D(i), time_since_periapsis/3600);
    fprintf('t- t_p at %d*R: %.7e days\n', D(i), time_since_periapsis/3600/24);
end;


% parabola
t = linspace(deg2rad(-140), deg2rad(140));
a = R_earth + 225;
x = -a*(t.^2) + (R_earth + 225);
y = 2*a*t;
plot(x,y);
hold on;

% center lines
yline(0);
xline(2*(R_earth + 225));

axis equal;
grid on;
hold off;
axis auto;
title("Parabolic Orbit of spacecraft around Earth")
xlabel("Distance [km]")
ylabel("Distance [km]")

% dominant earth-spacecraft
dominant = - (Gm_earth*75*R_earth) / (-75*R_earth)^3;
fprintf("dominant: %.7e km i^\n", dominant);

r_sc_moon = (75*R_earth) - (r_earth_moon);
direct = Gm_moon*(r_sc_moon/r_sc_moon^3);
indirect = Gm_moon*(-r_earth_moon/r_earth_moon^3);
net_perturbing = direct - indirect;
fprintf("direct: %.7e km i^\n", direct);
fprintf("indirect: %.7e km i^\n", indirect);
fprintf("net perturbing: %.7e km i^\n", net_perturbing);
