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

% Equations for Jupiter (Jupiter is now stationary in the rotating frame)
x_Jupiter = a_Jupiter*cos(0);
y_Jupiter = b_Jupiter*sin(0);

% The orbits of Asteroid A and Comet C relative to Jupiter for the entire orbit
results = zeros(0,5)
for i = 0:0.01:4
    t = i*Period_Jupiter;
    theta_comet = mean_motion_Comet*t;
    theta_asteroid = mean_motion_Asteroid*t;
    x_comet = a_comet * cos(theta_comet) - (a_comet*e);
    y_comet = b_comet * sin(theta_comet);
    TM = [cos(i*2*pi) sin(i*2*pi); -sin(i*2*pi) cos(i*2*pi)];
    comet_rotational = TM*[x_comet; y_comet];
    x_asteroid = a_asteroid * cos(theta_asteroid) - (a_asteroid*e);
    y_asteroid = b_asteroid * sin(theta_asteroid);
    asteroid_rotational = TM*[x_asteroid; y_asteroid];
    results(end+1,:) = [rad2deg(i*2*pi) comet_rotational(1) comet_rotational(2) asteroid_rotational(1) asteroid_rotational(2)];
end
t2 = array2table(results, 'VariableNames', {'Orbit', 'Comet X', 'Comet Y', 'Asteroid X', 'Asteroid Y'});
% disp(t2)

% Now, let's just get the positions when we're at 0, 0.25, 0.5 times Jupiter's period
results_select = zeros(0,7)
for i = 0:0.25:0.5
    t = i*Period_Jupiter;
    theta_comet = mean_motion_Comet*t;
    theta_asteroid = mean_motion_Asteroid*t;
    x_comet = a_comet * cos(theta_comet) - (a_comet*e);
    y_comet = b_comet * sin(theta_comet);
    TM = [cos(i*2*pi) sin(i*2*pi); -sin(i*2*pi) cos(i*2*pi)];
    comet_rotational = TM*[x_comet; y_comet];
    x_asteroid = a_asteroid * cos(theta_asteroid) - (a_asteroid*e);
    y_asteroid = b_asteroid * sin(theta_asteroid);
    asteroid_rotational = TM*[x_asteroid; y_asteroid];
    results_select(end+1,:) = [rad2deg(i*2*pi) rad2deg(theta_comet) rad2deg(theta_asteroid) comet_rotational(1) comet_rotational(2) asteroid_rotational(1) asteroid_rotational(2)];
end
t3 = array2table(results_select, 'VariableNames', {'Orbit', 'Theta_comet', 'Theta_asteroid', 'Comet_X', 'Comet_Y', 'Asteroid_X', 'Asteroid_Y'});
disp(t3)

% Plot just the orbits
fig1 = figure('Name', 'orbits only');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = scatter(x_Jupiter,y_Jupiter);
Comet_orbit = plot(t2, "Comet X", "Comet Y");
Asteroid_orbit = plot(t2, "Asteroid X", "Asteroid Y");
legend('Sun', 'Jupiter', 'Comet', 'Asteroid')
title('Orbits of Asteroid A and Comet C relative to the rotating frame')
xlim([-10e8 10e8])
ylim([-10e8 10e8])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot 0, 0.25, and 0.5.
fig1 = figure('Name', '0 and 0.25');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = scatter(x_Jupiter,y_Jupiter);
Comet_orbit = plot(t2, "Comet X", "Comet Y");
Asteroid_orbit = plot(t2, "Asteroid X", "Asteroid Y");
legend_0 = scatter(t3.Comet_X(1), t3.Comet_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
scatter(t3.Asteroid_X(1), t3.Asteroid_Y(1), 'black', 'square', 'filled', 'SizeData', 100)
legend_25 = scatter(t3.Comet_X(2), t3.Comet_Y(2), 'black', 'o', 'filled', 'SizeData', 100);
scatter(t3.Asteroid_X(2), t3.Asteroid_Y(2), 'black', 'o', 'filled', 'SizeData', 100)
legend_5 = scatter(t3.Comet_X(3), t3.Comet_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100);
scatter(t3.Asteroid_X(3), t3.Asteroid_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100)
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_0, legend_25, legend_5], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 't = 0P_J','t = 0.25P_J', 't = 0.5P_J'})
% legend('Sun', 'Jupiter', 'Comet', 'Asteroid')
title('Locations of Asteroid A and Comet C relative to the rotating frame')
xlim([-10e8 10e8])
ylim([-10e8 10e8])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')


% Try plotting 0.5 and 0.75 together
% fig1 = figure('Name', '0.5 and 0.75');
% scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200)
% hold on
% Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
% Comet_orbit = plot(x_comet,y_comet);
% Asteroid_orbit = plot(x_asteroid,y_asteroid);
% legend_5 = scatter(x_Jupiter_2, y_Jupiter_2, 'black', 'diamond', 'filled', 'SizeData', 100);
% scatter(x_comet_2, y_comet_2, 'black', 'diamond', 'filled', 'SizeData', 100)
% scatter(x_asteroid_2, y_asteroid_2, 'black', 'diamond', 'filled', 'SizeData', 100)
% legend_75 = scatter(x_Jupiter_3, y_Jupiter_3, 'black', 'o', 'filled', 'SizeData', 100);
% scatter(x_comet_3, y_comet_3, 'black', 'o', 'filled', 'SizeData', 100)
% scatter(x_asteroid_3, y_asteroid_3, 'black', 'o', 'filled', 'SizeData', 100)
% legend([Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_5, legend_75], {'Jupiter', 'Comet', 'Asteroid', 't = 0.5P_J','t = 0.75P_J'})
% title({'Locations of Jupiter, Asteroid A, and Comet C'; 'when t=0.5Period_{Jupiter} and t=0.75Period_{Jupiter}'})
% axis square
% hold off
% Jupiter_orbit.LineWidth = 2;
% Comet_orbit.LineWidth = 2;
% Asteroid_orbit.LineWidth = 2;
% fontsize(14, 'points')
% xlabel('Distance [km]')
% ylabel('Distance [km]')

% % Plot 1
% fig1 = figure('Name', '1');
% scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200)
% hold on
% Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
% Comet_orbit = plot(x_comet,y_comet);
% Asteroid_orbit = plot(x_asteroid,y_asteroid);
% for_legend = scatter(x_Jupiter_4, y_Jupiter_4, 'black', 'o', 'filled', 'SizeData', 100);
% scatter(x_comet_4, y_comet_4, 'black', 'o', 'filled', 'SizeData', 100)
% scatter(x_asteroid_4, y_asteroid_4, 'black', 'o', 'filled', 'SizeData', 100)
% legend([Jupiter_orbit, Comet_orbit, Asteroid_orbit, for_legend], {'Jupiter', 'Comet', 'Asteroid', 't = 1P_J'})
% title({'Locations of Jupiter, Asteroid A, and Comet C'; 'when t=1Period_{Jupiter}'})
% axis square
% hold off
% Jupiter_orbit.LineWidth = 2;
% Comet_orbit.LineWidth = 2;
% Asteroid_orbit.LineWidth = 2;
% fontsize(14, 'points')
% xlabel('Distance [km]')
% ylabel('Distance [km]')