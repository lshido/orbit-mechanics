R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
Gm_moon = 4902.8005821478;
mu = Gm_earth + Gm_moon;
a = 7.5*R_earth;
b = 3*R_earth;


theta = linspace(0,4*pi,200);
r = (a*cos(theta) + b*sin(theta))/R_earth;

fig1 = figure('Name', 'problem 1');
r_plot = plot(r);
hold on
fontsize(14, "points")
r_plot.LineWidth = 2;
hold off