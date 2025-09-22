% Givens
AU = 149597870.7; % km
a_Jupiter = 5.2*AU; % kma_Jupiter = 149597870.7 % km
mu_Sun = 132712440017.99; % km^3/s^2
mu_Jupiter = 126712767.8578; % km^3/s^2
mu_comet = 0; % assume massless
mu_asteroid = 0; % assume massless
e = 0.9;

% Determine the mean motion of Jupiter
mean_motion_Jupiter = sqrt((mu_Jupiter+mu_Sun)/a_Jupiter^3);
Period_Jupiter = 2*pi/mean_motion_Jupiter;
% Determine the mean motion of asteroid A.
mean_motion_Asteroid = (3/2)*mean_motion_Jupiter;
% Determine the mean motion of comet C.
mean_motion_Comet = (3/4)*mean_motion_Jupiter;

a_comet = -((mu_Sun+mu_comet)/mean_motion_Comet^2)^(1/3);
a_asteroid = ((mu_Sun+mu_asteroid)/mean_motion_Asteroid^2)^(1/3);

b_Jupiter = a_Jupiter; % a = b because circle
b_comet = a_comet*sqrt(1-e^2);
b_asteroid = a_asteroid*sqrt(1-e^2);

% Equations for the Sun
x_Sun = cos(0);
y_Sun = sin(0);

%=======================Plot the orbits=====================
% The orbits of Asteroid A and Comet C in the inertial frame
inertial_results = zeros(0,7);
for i = 0:0.01:2
    t = i*Period_Jupiter;
    % Equations for Jupiter
    x_Jupiter = a_Jupiter*cos(i*2*pi);
    y_Jupiter = b_Jupiter*sin(i*2*pi);
    % Equations for Comet
    theta_comet = mean_motion_Comet*t;
    x_comet = a_comet * cos(theta_comet) - (a_comet*e);
    y_comet = b_comet * sin(theta_comet);
    % Equations for Asteroid
    theta_asteroid = mean_motion_Asteroid*t;
    x_asteroid = a_asteroid * cos(theta_asteroid) - (a_asteroid*e);
    y_asteroid = b_asteroid * sin(theta_asteroid);
    inertial_results(end+1,:) = [rad2deg(i*2*pi) x_Jupiter y_Jupiter x_comet y_comet x_asteroid y_asteroid];
end
inertial_table = array2table(inertial_results, 'VariableNames', {'Orbit', 'Jupiter_X', 'Jupiter_Y', 'Comet_X', 'Comet_Y', 'Asteroid_X', 'Asteroid_Y'});

fig_anim = figure('Name', 'Animation!');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = plot(inertial_table, 'Jupiter_X', 'Jupiter_Y');
Comet_orbit = plot(inertial_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(inertial_table, "Asteroid_X", "Asteroid_Y");
Jupiter_marker = scatter(NaN,NaN,[], "black", "filled");
comet_marker = scatter(NaN,NaN,[], "black", "filled");
asteroid_marker = scatter(NaN,NaN,[], "black", "filled");
hold off
legend('Sun', 'Jupiter', 'Comet', 'Asteroid')
title({'Orbits of Jupiter, Asteroid A, and Comet C';['in the inertial frame with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

%=============ANIMATE THE POSITIONS!==========================

for i = 0:0.02:12
    t = i*Period_Jupiter;
    % Equations for Jupiter
    Jupiter_marker.XData = a_Jupiter*cos(i*2*pi);
    Jupiter_marker.YData = b_Jupiter*sin(i*2*pi);
    % Equations for Comet
    theta_comet = mean_motion_Comet*t;
    comet_marker.XData = a_comet * cos(theta_comet) - (a_comet*e);
    comet_marker.YData = b_comet * sin(theta_comet);
    % Equations for Asteroid
    theta_asteroid = mean_motion_Asteroid*t;
    asteroid_marker.XData = a_asteroid * cos(theta_asteroid) - (a_asteroid*e);
    asteroid_marker.YData = b_asteroid * sin(theta_asteroid);
    drawnow
    % exportgraphics(gca,"parabola.gif",Append=true)
end