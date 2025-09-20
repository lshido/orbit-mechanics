% Givens
AU = 149597870.7; % km
a_Jupiter = 5.2*AU; % kma_Jupiter = 149597870.7 % km
mu_Sun = 132712440017.99; % km^3/s^2
mu_Jupiter = 126712767.8578; % km^3/s^2
mu_comet = 0; % assume massless
mu_asteroid = 0; % assume massless

% Determine the mean motion of Jupiter
mean_motion_Jupiter = sqrt((mu_Jupiter+mu_Sun)/a_Jupiter^3);
Period_Jupiter = 2*pi/mean_motion_Jupiter;
% Determine the mean motion of asteroid A.
mean_motion_Asteroid = (3/2)*mean_motion_Jupiter;
% Determine the mean motion of comet C.
mean_motion_Comet = (3/4)*mean_motion_Jupiter;

a_comet = ((mu_Sun+mu_comet)/mean_motion_Comet^2)^(1/3);
a_asteroid = ((mu_Sun+mu_asteroid)/mean_motion_Asteroid^2)^(1/3);
fprintf('orbit radius of comet: %.9e [km]\n', a_comet);
fprintf('orbit radius of asteroid: %.9e [km]\n', a_asteroid);

% Calculate the theta value for the asteroid and comet based on Jupiter's position
% Step 1: Initialize an array
% Step 2: For every 0.25 increment in Jupiter's period
% Step 3: Calculate the time
% Step 4: Calculate the theta angle for the object
results_theta = zeros(0,4);
for i = 0:0.25:1
    t = i*Period_Jupiter;
    theta_comet = mean_motion_Comet*t;
    theta_asteroid = mean_motion_Asteroid*t;
    % maybe plot here...?
    results_theta(end+1,:) = [i rad2deg(i*2*pi) mod(rad2deg(theta_comet), 360) mod(rad2deg(theta_asteroid), 360)];
end

t1 = array2table(results_theta, 'VariableNames', {'Fraction of Jupiter Period', 'Theta of Jupiter [deg]', 'Theta of Comet [deg]', 'Theta of Asteroid [deg]'});
disp(t1)

e = 0; % circle
b_Jupiter = a_Jupiter; % a = b because circle
b_comet = a_comet; % a = b because circle
b_asteroid = a_asteroid; % a = b because circle

% Equations for the Sun
x_Sun = cos(0);
y_Sun = sin(0);

% Equations for the orbits
t = linspace(0, 2*pi);
x_Jupiter = a_Jupiter * cos(t) - (a_Jupiter*e);
y_Jupiter = b_Jupiter * sin(t);

x_comet = a_comet * cos(t) - (a_comet*e);
y_comet = b_comet * sin(t);

x_asteroid = a_asteroid * cos(t) - (a_asteroid*e);
y_asteroid = b_asteroid * sin(t);

% Equations for the first position
t = 0;
x_Jupiter_0 = a_Jupiter * cos(t) - (a_Jupiter*e);
y_Jupiter_0 = b_Jupiter * sin(t);

x_comet_0 = a_comet * cos(t) - (a_comet*e);
y_comet_0 = b_comet * sin(t);

x_asteroid_0 = a_asteroid * cos(t) - (a_asteroid*e);
y_asteroid_0 = b_asteroid * sin(t);

% Equations for the second position
t_J = deg2rad(90); 
t_C = deg2rad(67.5);
t_A = deg2rad(135);
x_Jupiter_1 = a_Jupiter * cos(t_J) - (a_Jupiter*e);
y_Jupiter_1 = b_Jupiter * sin(t_J);

x_comet_1 = a_comet * cos(t_C) - (a_comet*e);
y_comet_1 = b_comet * sin(t_C);

x_asteroid_1 = a_asteroid * cos(t_A) - (a_asteroid*e);
y_asteroid_1 = b_asteroid * sin(t_A);

% Equations for the third position
t_J = deg2rad(180); 
t_C = deg2rad(135);
t_A = deg2rad(270);
x_Jupiter_2 = a_Jupiter * cos(t_J) - (a_Jupiter*e);
y_Jupiter_2 = b_Jupiter * sin(t_J);

x_comet_2 = a_comet * cos(t_C) - (a_comet*e);
y_comet_2 = b_comet * sin(t_C);

x_asteroid_2 = a_asteroid * cos(t_A) - (a_asteroid*e);
y_asteroid_2 = b_asteroid * sin(t_A);

% Equations for the fourth position
t_J = deg2rad(270); 
t_C = deg2rad(202.5);
t_A = deg2rad(45);
x_Jupiter_3 = a_Jupiter * cos(t_J) - (a_Jupiter*e);
y_Jupiter_3 = b_Jupiter * sin(t_J);

x_comet_3 = a_comet * cos(t_C) - (a_comet*e);
y_comet_3 = b_comet * sin(t_C);

x_asteroid_3 = a_asteroid * cos(t_A) - (a_asteroid*e);
y_asteroid_3 = b_asteroid * sin(t_A);

% Equations for the fifth position
t_J = deg2rad(360); 
t_C = deg2rad(270);
t_A = deg2rad(180);
x_Jupiter_4 = a_Jupiter * cos(t_J) - (a_Jupiter*e);
y_Jupiter_4 = b_Jupiter * sin(t_J);

x_comet_4 = a_comet * cos(t_C) - (a_comet*e);
y_comet_4 = b_comet * sin(t_C);

x_asteroid_4 = a_asteroid * cos(t_A) - (a_asteroid*e);
y_asteroid_4 = b_asteroid * sin(t_A);

% Plot 0 and 0.25
fig1 = figure('Name', '0 and 0.25');
scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200)
hold on
Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
Comet_orbit = plot(x_comet,y_comet);
Asteroid_orbit = plot(x_asteroid,y_asteroid);
legend_0 = scatter(x_Jupiter_0, y_Jupiter_0, 'black', 'square', 'filled', 'SizeData', 100);
scatter(x_comet_0, y_comet_0, 'black', 'square', 'filled', 'SizeData', 100)
scatter(x_asteroid_0, y_asteroid_0, 'black', 'square', 'filled', 'SizeData', 100)
legend_25 = scatter(x_Jupiter_1, y_Jupiter_1, 'black', 'o', 'filled', 'SizeData', 100);
scatter(x_comet_1, y_comet_1, 'black', 'o', 'filled', 'SizeData', 100)
scatter(x_asteroid_1, y_asteroid_1, 'black', 'o', 'filled', 'SizeData', 100)
legend([Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_0, legend_25], {'Jupiter', 'Comet', 'Asteroid', 't = 0P_J','t = 0.25P_J'})
title({'Locations of Jupiter, Asteroid A, and Comet C'; 'when t=0Period_{Jupiter} and t=0.25Period_{Jupiter}'})
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')


% Try plotting 0.5 and 0.75 together
fig1 = figure('Name', '0.5 and 0.75');
scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200)
hold on
Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
Comet_orbit = plot(x_comet,y_comet);
Asteroid_orbit = plot(x_asteroid,y_asteroid);
legend_5 = scatter(x_Jupiter_2, y_Jupiter_2, 'black', 'diamond', 'filled', 'SizeData', 100);
scatter(x_comet_2, y_comet_2, 'black', 'diamond', 'filled', 'SizeData', 100)
scatter(x_asteroid_2, y_asteroid_2, 'black', 'diamond', 'filled', 'SizeData', 100)
legend_75 = scatter(x_Jupiter_3, y_Jupiter_3, 'black', 'o', 'filled', 'SizeData', 100);
scatter(x_comet_3, y_comet_3, 'black', 'o', 'filled', 'SizeData', 100)
scatter(x_asteroid_3, y_asteroid_3, 'black', 'o', 'filled', 'SizeData', 100)
legend([Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_5, legend_75], {'Jupiter', 'Comet', 'Asteroid', 't = 0.5P_J','t = 0.75P_J'})
title({'Locations of Jupiter, Asteroid A, and Comet C'; 'when t=0.5Period_{Jupiter} and t=0.75Period_{Jupiter}'})
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot 1
fig1 = figure('Name', '1');
scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200)
hold on
Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
Comet_orbit = plot(x_comet,y_comet);
Asteroid_orbit = plot(x_asteroid,y_asteroid);
for_legend = scatter(x_Jupiter_4, y_Jupiter_4, 'black', 'o', 'filled', 'SizeData', 100);
scatter(x_comet_4, y_comet_4, 'black', 'o', 'filled', 'SizeData', 100)
scatter(x_asteroid_4, y_asteroid_4, 'black', 'o', 'filled', 'SizeData', 100)
legend([Jupiter_orbit, Comet_orbit, Asteroid_orbit, for_legend], {'Jupiter', 'Comet', 'Asteroid', 't = 1P_J'})
title({'Locations of Jupiter, Asteroid A, and Comet C'; 'when t=1Period_{Jupiter}'})
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')