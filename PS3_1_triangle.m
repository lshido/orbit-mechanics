r_sun_jupiter = [ +778279959 0 ];
theta = 60;  % same for all angles

% Problem 1b(i)
mag_r_sun_asteroid = norm(r_sun_jupiter) / (2 * cosd(theta));
fprintf('magnitude of r_sun_asteroid: %.2e\n', mag_r_sun_asteroid);

i_r_sun_asteroid = mag_r_sun_asteroid * cosd(theta);
j_r_sun_asteroid = mag_r_sun_asteroid * sind(theta);
fprintf('r_sun_asteroid: %.2ei %.2ej\n', i_r_sun_asteroid, j_r_sun_asteroid);

mag_r_asteroid_jupiter = norm(r_sun_jupiter) / (2 * cosd(theta));

i_r_asteroid_jupiter = mag_r_asteroid_jupiter * cosd(theta);
j_r_asteroid_jupiter = mag_r_asteroid_jupiter * -sind(theta);
fprintf('r_asteroid_jupiter: %.2ei %.2ej\n', i_r_asteroid_jupiter, j_r_asteroid_jupiter);