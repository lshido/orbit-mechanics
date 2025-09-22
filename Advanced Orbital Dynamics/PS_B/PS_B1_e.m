% Givens
AU = 149597870.7; % km
a_Jupiter = 5.2*AU; % kma_Jupiter = 149597870.7 % km
mu_Sun = 132712440017.99; % km^3/s^2
mu_Jupiter = 126712767.8578; % km^3/s^2
mu_comet = 0; % assume massless
mu_asteroid = 0; % assume massless
e = 0.55;

% Determine the mean motion of Jupiter
mean_motion_Jupiter = sqrt((mu_Jupiter+mu_Sun)/a_Jupiter^3);
Period_Jupiter = 2*pi/mean_motion_Jupiter;
% Determine the mean motion of asteroid A.
mean_motion_Asteroid = (3/2)*mean_motion_Jupiter;
% Determine the mean motion of comet C.
mean_motion_Comet = (3/4)*mean_motion_Jupiter;

a_comet = -((mu_Sun+mu_comet)/mean_motion_Comet^2)^(1/3);
a_asteroid = ((mu_Sun+mu_asteroid)/mean_motion_Asteroid^2)^(1/3);

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

b_Jupiter = a_Jupiter; % a = b because circle
b_comet = a_comet*sqrt(1-e^2);
b_asteroid = a_asteroid*sqrt(1-e^2);

% Equations for the Sun
x_Sun = cos(0);
y_Sun = sin(0);

%================PLOT THE INERTIAL FRAME==============================

% The orbits of Asteroid A and Comet C in the inertial frame
inertial_results = zeros(0,7);
for i = 0:0.01:4
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
% disp(inertial_table)

% Now, let's just get the positions when we're at 0, 0.25, 0.5 times Jupiter's period
inertial_results_select = zeros(0,7);
for i = 0:0.25:1
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
    inertial_results_select(end+1,:) = [rad2deg(i*2*pi) x_Jupiter y_Jupiter x_comet y_comet x_asteroid y_asteroid];
end
inertial_select = array2table(inertial_results_select, 'VariableNames', {'Orbit', 'Jupiter_X', 'Jupiter_Y', 'Comet_X', 'Comet_Y', 'Asteroid_X', 'Asteroid_Y'});
disp(inertial_select)

% Plot just the inertial orbits
fig1 = figure('Name', 'inertial orbits');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = plot(inertial_table, 'Jupiter_X', 'Jupiter_Y');
Comet_orbit = plot(inertial_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(inertial_table, "Asteroid_X", "Asteroid_Y");
legend('Sun', 'Jupiter', 'Comet', 'Asteroid')
title({'Orbits of Jupiter, Asteroid A, and Comet C';['in the inertial frame with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot the 0, 0.25, 0.5 positions in the inertial frame
fig1 = figure('Name', 'inertial Positions');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = plot(inertial_table, 'Jupiter_X', 'Jupiter_Y');
Comet_orbit = plot(inertial_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(inertial_table, "Asteroid_X", "Asteroid_Y");
legend_0 = scatter(inertial_select.Jupiter_X(1), inertial_select.Jupiter_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
scatter(inertial_select.Comet_X(1), inertial_select.Comet_Y(1), 'black', 'square', 'filled', 'SizeData', 100)
scatter(inertial_select.Asteroid_X(1), inertial_select.Asteroid_Y(1), 'black', 'square', 'filled', 'SizeData', 100)
legend_25 = scatter(inertial_select.Jupiter_X(2), inertial_select.Jupiter_Y(2), 'black', 'o', 'filled', 'SizeData', 100);
scatter(inertial_select.Comet_X(2), inertial_select.Comet_Y(2), 'black', 'o', 'filled', 'SizeData', 100)
scatter(inertial_select.Asteroid_X(2), inertial_select.Asteroid_Y(2), 'black', 'o', 'filled', 'SizeData', 100)
legend_5 = scatter(inertial_select.Jupiter_X(3), inertial_select.Jupiter_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100);
scatter(inertial_select.Comet_X(3), inertial_select.Comet_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100)
scatter(inertial_select.Asteroid_X(3), inertial_select.Asteroid_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100)
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_0, legend_25, legend_5], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 't = 0P_J','t = 0.25P_J', 't = 0.5P_J'})
title({'Locations of Jupiter, Asteroid A, and Comet C';['in the inertial frame at various times with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot 0.75 and 1 positions in inertial frame
fig1 = figure('Name', 'inertial Positions');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = plot(inertial_table, 'Jupiter_X', 'Jupiter_Y');
Comet_orbit = plot(inertial_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(inertial_table, "Asteroid_X", "Asteroid_Y");
legend_75 = scatter(inertial_select.Jupiter_X(4), inertial_select.Jupiter_Y(4), 'black', 'square', 'filled', 'SizeData', 100);
scatter(inertial_select.Comet_X(4), inertial_select.Comet_Y(4), 'black', 'square', 'filled', 'SizeData', 100)
scatter(inertial_select.Asteroid_X(4), inertial_select.Asteroid_Y(4), 'black', 'square', 'filled', 'SizeData', 100)
legend_1 = scatter(inertial_select.Jupiter_X(5), inertial_select.Jupiter_Y(5), 'black', 'o', 'filled', 'SizeData', 100);
scatter(inertial_select.Comet_X(5), inertial_select.Comet_Y(5), 'black', 'o', 'filled', 'SizeData', 100)
scatter(inertial_select.Asteroid_X(5), inertial_select.Asteroid_Y(5), 'black', 'o', 'filled', 'SizeData', 100)
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, legend_75, legend_1], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 't = 0.75P_J','t = 1P_J'})
title({'Locations of Jupiter, Asteroid A, and Comet C';['in the inertial frame at various times with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

%=================PLOT THE ROTATIONAL FRAME===============================
% Equations for Jupiter (Jupiter is now stationary in the rotating frame)
x_Jupiter = a_Jupiter*cos(0);
y_Jupiter = b_Jupiter*sin(0);

% The orbits of Asteroid A and Comet C relative to Jupiter for the entire orbit
rotational_results = zeros(0,5)
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
    rotational_results(end+1,:) = [rad2deg(i*2*pi) comet_rotational(1) comet_rotational(2) asteroid_rotational(1) asteroid_rotational(2)];
end
rotational_table = array2table(rotational_results, 'VariableNames', {'Orbit', 'Comet_X', 'Comet_Y', 'Asteroid_X', 'Asteroid_Y'});
% disp(rotational_table)

% Now, let's just get the positions when we're at 0, 0.25, 0.5 times Jupiter's period
rotational_results_select = zeros(0,5)
for i = 0:0.25:1
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
    rotational_results_select(end+1,:) = [rad2deg(i*2*pi) comet_rotational(1) comet_rotational(2) asteroid_rotational(1) asteroid_rotational(2)];
end
rotational_select = array2table(rotational_results_select, 'VariableNames', {'Orbit', 'Comet_X', 'Comet_Y', 'Asteroid_X', 'Asteroid_Y'});
disp(rotational_select)

% Get the periapsis x,y values in the rotating view. 
peri_results = zeros(0,7)
for angle = [0 360 720 1080 1440]
    i_c = deg2rad(angle)/(mean_motion_Comet*Period_Jupiter);
    TM_c = [cos(i_c*2*pi) sin(i_c*2*pi); -sin(i_c*2*pi) cos(i_c*2*pi)];
    x_comet = a_comet * cos(deg2rad(angle)) - (a_comet*e);
    y_comet = b_comet * sin(deg2rad(angle));
    comet_rotational = TM_c*[x_comet; y_comet];
    i_a = deg2rad(angle)/(mean_motion_Asteroid*Period_Jupiter);
    TM_a = [cos(i_a*2*pi) sin(i_a*2*pi); -sin(i_a*2*pi) cos(i_a*2*pi)];
    x_asteroid = a_asteroid * cos(deg2rad(angle)) - (a_asteroid*e);
    y_asteroid = b_asteroid * sin(deg2rad(angle));
    asteroid_rotational = TM_a*[x_asteroid; y_asteroid];
    peri_results(end+1,:) = [angle i_c comet_rotational(1) comet_rotational(2) i_a asteroid_rotational(1) asteroid_rotational(2)];
end
peri_table = array2table(peri_results, 'VariableNames', {'Periapsis', 'i_c', 'Comet_X', 'Comet_Y', 'i_a', 'Asteroid_X', 'Asteroid_Y'});
disp(peri_table)

% Get the apoapsis x,y values in the rotating view. 
apo_results = zeros(0,7)
for angle = [180 540 900 1260 1620]
    i_c = deg2rad(angle)/(mean_motion_Comet*Period_Jupiter);
    TM_c = [cos(i_c*2*pi) sin(i_c*2*pi); -sin(i_c*2*pi) cos(i_c*2*pi)];
    x_comet = a_comet * cos(deg2rad(angle)) - (a_comet*e);
    y_comet = b_comet * sin(deg2rad(angle));
    comet_rotational = TM_c*[x_comet; y_comet];
    i_a = deg2rad(angle)/(mean_motion_Asteroid*Period_Jupiter);
    TM_a = [cos(i_a*2*pi) sin(i_a*2*pi); -sin(i_a*2*pi) cos(i_a*2*pi)];
    x_asteroid = a_asteroid * cos(deg2rad(angle)) - (a_asteroid*e);
    y_asteroid = b_asteroid * sin(deg2rad(angle));
    asteroid_rotational = TM_a*[x_asteroid; y_asteroid];
    apo_results(end+1,:) = [angle rad2deg(i_c*2*pi) comet_rotational(1) comet_rotational(2) rad2deg(i_a*2*pi) asteroid_rotational(1) asteroid_rotational(2)];
end
apo_table = array2table(apo_results, 'VariableNames', {'Apoapsis', 'Angle of rotational frame x axis at comet apo', 'Comet_X', 'Comet_Y', 'Angle of rotational frame x axis at asteroid apo', 'Asteroid_X', 'Asteroid_Y'});
disp(apo_table)

% Plot just the rotational orbits
fig1 = figure('Name', 'rotational orbits');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = scatter(x_Jupiter,y_Jupiter);
Comet_orbit = plot(rotational_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(rotational_table, "Asteroid_X", "Asteroid_Y");
% Comets peri and apo
comet_Peri1_legend = scatter(peri_table.Comet_X(1),peri_table.Comet_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
comet_Peri2_legend = scatter(peri_table.Comet_X(2),peri_table.Comet_Y(2), 'black', 'square', 'filled', 'SizeData', 100);
comet_Peri3_legend = scatter(peri_table.Comet_X(3),peri_table.Comet_Y(3), 'black', 'square', 'filled', 'SizeData', 100);
comet_Apo1_legend = scatter(apo_table.Comet_X(1),apo_table.Comet_Y(1), 'black', 'filled', 'SizeData', 100);
comet_Apo2_legend = scatter(apo_table.Comet_X(2),apo_table.Comet_Y(2), 'black', 'filled', 'SizeData', 100);
comet_Apo3_legend = scatter(apo_table.Comet_X(3),apo_table.Comet_Y(3), 'black', 'filled', 'SizeData', 100);
% Asteroid peri and apo
asteroid_Peri1 = scatter(peri_table.Asteroid_X(1),peri_table.Asteroid_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
asteroid_Peri2 = scatter(peri_table.Asteroid_X(2),peri_table.Asteroid_Y(2), 'black', 'square', 'filled', 'SizeData', 100);
asteroid_Peri3 = scatter(peri_table.Asteroid_X(3),peri_table.Asteroid_Y(3), 'black', 'square', 'filled', 'SizeData', 100);
asteroid_Apo1 = scatter(apo_table.Asteroid_X(1),apo_table.Asteroid_Y(1), 'black', 'filled', 'SizeData', 100);
asteroid_Apo2 = scatter(apo_table.Asteroid_X(2),apo_table.Asteroid_Y(2), 'black', 'filled', 'SizeData', 100);
asteroid_Apo3 = scatter(apo_table.Asteroid_X(3),apo_table.Asteroid_Y(3), 'black', 'filled', 'SizeData', 100);
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, comet_Peri1_legend, comet_Apo1_legend], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 'Periapsis', 'Apoapsis'})
title({'Orbits of Asteroid A and Comet C';['relative to the rotating frame with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot 0, 0.25, and 0.5.
fig1 = figure('Name', 'rotational positions');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = scatter(x_Jupiter,y_Jupiter);
Comet_orbit = plot(rotational_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(rotational_table, "Asteroid_X", "Asteroid_Y");
comet_0_legend = scatter(rotational_select.Comet_X(1), rotational_select.Comet_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
asteroid_0 = scatter(rotational_select.Asteroid_X(1), rotational_select.Asteroid_Y(1), 'black', 'square', 'filled', 'SizeData', 100);
comet_25_legend = scatter(rotational_select.Comet_X(2), rotational_select.Comet_Y(2), 'black', 'o', 'filled', 'SizeData', 100);
asteroid_25 = scatter(rotational_select.Asteroid_X(2), rotational_select.Asteroid_Y(2), 'black', 'o', 'filled', 'SizeData', 100);
comet_5_legend = scatter(rotational_select.Comet_X(3), rotational_select.Comet_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100);
asteroid_5 = scatter(rotational_select.Asteroid_X(3), rotational_select.Asteroid_Y(3), 'black', 'diamond', 'filled', 'SizeData', 100);
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, comet_0_legend, comet_25_legend, comet_5_legend], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 't = 0P_J','t = 0.25P_J', 't = 0.5P_J'})
title({'Orbits of Asteroid A and Comet C';['relative to the rotating frame with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')

% Plot 0.75 and 1 in rotational frame.
fig1 = figure('Name', 'rotational positions');
sun = scatter(x_Sun, y_Sun, 'red', 'filled', 'SizeData', 200);
hold on
Jupiter_orbit = scatter(x_Jupiter,y_Jupiter);
Comet_orbit = plot(rotational_table, "Comet_X", "Comet_Y");
Asteroid_orbit = plot(rotational_table, "Asteroid_X", "Asteroid_Y");
comet_75_legend = scatter(rotational_select.Comet_X(4), rotational_select.Comet_Y(4), 'black', 'square', 'filled', 'SizeData', 100);
asteroid_75 = scatter(rotational_select.Asteroid_X(4), rotational_select.Asteroid_Y(4), 'black', 'square', 'filled', 'SizeData', 100);
comet_1_legend = scatter(rotational_select.Comet_X(5), rotational_select.Comet_Y(5), 'black', 'o', 'filled', 'SizeData', 100);
asteroid_1 = scatter(rotational_select.Asteroid_X(5), rotational_select.Asteroid_Y(5), 'black', 'o', 'filled', 'SizeData', 100);
legend([sun, Jupiter_orbit, Comet_orbit, Asteroid_orbit, comet_75_legend, comet_1_legend], {'Sun', 'Jupiter', 'Comet', 'Asteroid', 't = 0.75P_J','t = 1P_J'})
title({'Locations of Asteroid A and Comet C';['relative to the rotating frame with e=',num2str(e)]})
xlim([-1.5e9 1.5e9])
ylim([-1.5e9 1.5e9])
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, 'points')
xlabel('Distance [km]')
ylabel('Distance [km]')