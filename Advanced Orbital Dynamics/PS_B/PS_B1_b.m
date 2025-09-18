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
fprintf("orbit radius of comet: %.9e [km]\n", a_comet);
fprintf("orbit radius of asteroid: %.9e [km]\n", a_asteroid);

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

% Equations for the orbits
t = linspace(0, 2*pi);
x_Jupiter = a_Jupiter * cos(t) - (a_Jupiter*e);
y_Jupiter = b_Jupiter * sin(t);

x_comet = a_comet * cos(t) - (a_comet*e);
y_comet = b_comet * sin(t);

x_asteroid = a_asteroid * cos(t) - (a_asteroid*e);
y_asteroid = b_asteroid * sin(t);

% Equations for the first position
% TBD

% Plot it
Jupiter_orbit = plot(x_Jupiter,y_Jupiter);
hold on
Comet_orbit = plot(x_comet,y_comet);
Asteroid_orbit = plot(x_asteroid,y_asteroid);
legend("Jupiter", "Comet", "Asteroid")
axis square
hold off
Jupiter_orbit.LineWidth = 2;
Comet_orbit.LineWidth = 2;
Asteroid_orbit.LineWidth = 2;
fontsize(14, "points")
xlabel("Distance [km]")
ylabel("Distance [km]")