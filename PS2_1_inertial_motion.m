% known position vectors
r_moon_sc = +140000;
r_sun_earth = -149597898;
r_earth_moon = +384400;
r_sun_jupiter = +778279959;

% known Gm values
Gm_earth = 398600.4415;
Gm_sun = 132712440017.99;
Gm_moon = 4902.8005821478;
Gm_jupiter = 126712767.8578;


% Acceleration due to earth
r_earth_sc = r_earth_moon + r_moon_sc;
accel_earth = -Gm_earth * r_earth_sc / (abs(r_earth_sc)^3);
fprintf('accel_earth = %.5e\n', accel_earth);

% Acceleration due to moon
accel_moon = -Gm_moon * r_moon_sc / (abs(r_moon_sc)^3);
fprintf('accel_moon = %.5e\n', accel_moon);

% Acceleration due to sun
r_sun_sc = r_sun_earth + r_earth_moon + r_moon_sc;
accel_sun = -Gm_sun * r_sun_sc / (abs(r_sun_sc)^3);
fprintf('accel_sun = %.5e\n', accel_sun);

% Acceleration due to jupiter
r_jupiter_sc = r_sun_earth + r_earth_moon + r_moon_sc - r_sun_jupiter;
accel_jupiter = -Gm_jupiter * r_jupiter_sc / (abs(r_jupiter_sc)^3);
fprintf('accel_jupiter = %.5e\n', accel_jupiter);

net_accel = accel_earth + accel_moon + accel_sun + accel_jupiter;
fprintf('net_accel = %.5e\n', net_accel);

