R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
Gm_moon = 4902.8005821478;
mu =1e-4;
% mu = Gm_earth + Gm_moon;
a = 7.5*R_earth;
b = 3*R_earth;

% theta= linspace(0,360,360);
% x = a*cos(deg2rad(theta));
% y = b*sin(deg2rad(theta));
% x1 = 10*a*cos(deg2rad(theta));
% y1 = b*sin(deg2rad(theta));
% x2 = a*cos(deg2rad(theta));
% y2 = 10*b*sin(deg2rad(theta));
% x3 = a*cos(deg2rad(theta));
% y3 = a*sin(deg2rad(theta));

t = linspace(0,3600,3600);
x = a*cos(sqrt(mu)*t);
y = b*sin(sqrt(mu)*t);
x1 = 10*a*cos(sqrt(mu)*t);
y1 = b*sin(sqrt(mu)*t);
x2 = a*cos(sqrt(mu)*t);
y2 = 10*b*sin(sqrt(mu)*t);
x3 = a*cos(sqrt(mu)*t);
y3 = a*sin(sqrt(mu)*t);

r = sqrt(x.^2 + y.^2)/R_earth;
r1 = sqrt(x1.^2 + y1.^2)/R_earth;
r2 = sqrt(x2.^2 + y2.^2)/R_earth;
r3 = sqrt(x3.^2 + y3.^2)/R_earth;

% r5 = (x + y)/R_earth;
% r6 = (x1 + y1)/R_earth;
% r7 = (x2 + y2)/R_earth;
% r8 = (x3 + y3)/R_earth;

s = 2000;
beta = linspace(0, 2*pi);
my_dot_x = s*cos(beta);
my_dot_y = s*sin(beta);

fig1 = figure('Name', 'problem 1-orbits');
x_plot = plot(x,y);
axis square
hold on
x1_plot = plot(x1,y1);
x2_plot = plot(x2,y2);
x3_plot = plot(x3,y3);
my_dot = plot(my_dot_x,my_dot_y);

fontsize(14, "points")
x_plot.LineWidth = 2;
x1_plot.LineWidth = 2;
x2_plot.LineWidth = 2;
x3_plot.LineWidth = 2;
my_dot.LineWidth = 2;
legend("a=7.5R_{earth},b=3*R_{earth}","a=75R_{earth},b=3*R_{earth}", "a=7.5R_{earth},b=30R_{earth}", "a=7.5R_{earth},b=7.5R_{earth}")
ylabel("Y-Position in inertial frame [km]")
xlabel("X-Position in inertial frame [km]")
title("Orbits represented by r vector (Lillian Shido, PSA\_1\_e.m)")

xlim([-5e5 5e5])
ylim([-5e5 5e5])
hold off

fig2 = figure('Name','problem1-magnitude of r');
r_plot = plot(r);
hold on
r1_plot = plot(r1);
r2_plot = plot(r2);
r3_plot = plot(r3);
% r4_plot = plot(r4);
fontsize(14, "points")

r_plot.LineWidth = 2;
r1_plot.LineWidth = 2;
r2_plot.LineWidth = 2;
r3_plot.LineWidth = 2;
% r4_plot.LineWidth = 2;
ylabel("Magnitude of r in R_{earth}")
% xlabel("Values for \theta in degrees")
xlabel("Time in seconds")
title("Relative motion (Lillian Shido, PSA\_1\_e.m)")
legend("a=7.5R_{earth},b=3*R_{earth}","a=75R_{earth},b=3*R_{earth}", "a=7.5R_{earth},b=30R_{earth}", "a=7.5R_{earth},b=7.5R_{earth}")
% xlim([0 360])
xlim([0 3600])
hold off

% Various mus
mu1 =1e-6;
xm1 = a*cos(sqrt(mu1)*t);
ym1 = b*sin(sqrt(mu1)*t);
rm1 = sqrt(xm1.^2 + ym1.^2)/R_earth;

mu2 =1e-5;
xm2 = a*cos(sqrt(mu2)*t);
ym2 = b*sin(sqrt(mu2)*t);
rm2 = sqrt(xm2.^2 + ym2.^2)/R_earth;

mu3 =1e-4;
xm3 = a*cos(sqrt(mu3)*t);
ym3 = b*sin(sqrt(mu3)*t);
rm3 = sqrt(xm3.^2 + ym3.^2)/R_earth;

fig3 = figure('Name','problem1-r motion');
rm1_plot = plot(rm1);
hold on
rm2_plot = plot(rm2);
rm3_plot = plot(rm3);
% r4_plot = plot(r4);
fontsize(14, "points")

rm1_plot.LineWidth = 2;
rm2_plot.LineWidth = 2;
rm3_plot.LineWidth = 2;
% r4_plot.LineWidth = 2;
ylabel("Magnitude of r in R_{earth}")
% xlabel("Values for \theta in degrees")
xlabel("Time in seconds")
title("Relative motion with different \mu values (Lillian Shido, PSA\_1\_e.m)")
legend("\mu=1e-6","\mu=1e-5", "\mu=1e-4")
% xlim([0 360])
xlim([0 3600])
hold off