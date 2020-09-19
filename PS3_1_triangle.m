r_sun_jupiter = [ +778279959 0 ];
theta = 60;  % same for all angles
Gm_asteroid = 75;
Gm_sun = 132712440017.99;
Gm_jupiter = 126712767.8578;

% Problem 1b(i)
mag_r_sun_asteroid = norm(r_sun_jupiter);
fprintf('magnitude of r_sun_asteroid: %.2e\n', mag_r_sun_asteroid);

i_r_sun_asteroid = mag_r_sun_asteroid * cosd(theta);
j_r_sun_asteroid = mag_r_sun_asteroid * sind(theta);
r_sun_asteroid = [i_r_sun_asteroid j_r_sun_asteroid];
fprintf('r_sun_asteroid: %.2ei %.2ej\n', r_sun_asteroid(1), r_sun_asteroid(2));

mag_r_asteroid_jupiter = norm(r_sun_jupiter);
fprintf('magnitude of r_asteroid_jupiter: %.2e\n', mag_r_asteroid_jupiter);

i_r_asteroid_jupiter = mag_r_asteroid_jupiter * cosd(theta);
j_r_asteroid_jupiter = mag_r_asteroid_jupiter * -sind(theta);
r_asteroid_jupiter = [i_r_asteroid_jupiter j_r_asteroid_jupiter];
fprintf('r_asteroid_jupiter: %.2ei %.2ej\n', r_asteroid_jupiter(1), r_asteroid_jupiter(2));


% calculate dominant term for asteroid relative to sun
dominant_sun_asteroid = ((-Gm_asteroid - Gm_sun) * r_sun_asteroid)/((norm(r_sun_asteroid))^3);
fprintf('dominant_sun_asteroid = %.5ei %.5ej km/s^2\n', dominant_sun_asteroid(1), dominant_sun_asteroid(2));
fprintf('mag_dominant_sun_asteroid = %.5e km/s^2\n', norm(dominant_sun_asteroid));

% perturbing due to jupiter
direct_jupiter = Gm_jupiter*r_asteroid_jupiter/(norm(r_asteroid_jupiter))^3;
indirect_jupiter = Gm_jupiter*r_sun_jupiter/(norm(r_sun_jupiter))^3;
net_pert_jupiter = direct_jupiter - indirect_jupiter;
fprintf('direct_jupiter = %.5ei %.5ej km/s^2\n', direct_jupiter(1), direct_jupiter(2));
fprintf('mag_direct_jupiter = %.5e km/s^2\n', norm(direct_jupiter));
fprintf('indirect_jupiter = %.5ei %.5ej km/s^2\n', indirect_jupiter(1), indirect_jupiter(2));
fprintf('mag_indirect_jupiter = %.5e km/s^2\n', norm(indirect_jupiter));
fprintf('net_pert_jupiter = %.5ei %.5ej km/s^2\n', net_pert_jupiter(1), net_pert_jupiter(2));
fprintf('mag_net_pert_jupiter = %.5e km/s^2\n', norm(net_pert_jupiter));

total_net_sun_asteroid = dominant_sun_asteroid + net_pert_jupiter;
fprintf('total net sun asteroid = %.5ei %.5ej km/s^2\n', total_net_sun_asteroid(1), total_net_sun_asteroid(2));
fprintf('mag total net sun asteroid = %.5e km/s^2\n\n', norm(total_net_sun_asteroid));



% asteroid relative to jupiter
mag_r_asteroid_sun = norm(r_sun_jupiter);
fprintf('magnitude of r_sun_asteroid: %.2e\n', mag_r_sun_asteroid);

i_r_asteroid_sun = -mag_r_asteroid_sun * cosd(theta);
j_r_asteroid_sun = -mag_r_asteroid_sun * sind(theta);
r_asteroid_sun = [i_r_asteroid_sun j_r_asteroid_sun];
fprintf('r_asteroid_sun: %.2ei %.2ej\n', r_asteroid_sun(1), r_asteroid_sun(2));

mag_r_jupiter_asteroid = norm(r_sun_jupiter);
fprintf('magnitude of r_jupiter_asteroid: %.2e\n', mag_r_jupiter_asteroid);

i_r_jupiter_asteroid = -mag_r_jupiter_asteroid * cosd(theta);
j_r_jupiter_asteroid = mag_r_jupiter_asteroid * sind(theta);
r_jupiter_asteroid = [i_r_jupiter_asteroid j_r_jupiter_asteroid];
fprintf('r_jupiter_asteroid: %.2ei %.2ej\n', r_jupiter_asteroid(1), r_jupiter_asteroid(2));


r_jupiter_sun = [-norm(r_sun_jupiter) 0];
fprintf('r_jupiter_sun: %.2ei %.2ej\n', r_jupiter_sun(1), r_jupiter_sun(2));

% calculate dominant term for asteroid relative to jupiter
dominant_jupiter_asteroid = ((-Gm_asteroid - Gm_jupiter) * r_jupiter_asteroid)/((norm(r_jupiter_asteroid))^3);
fprintf('dominant_jupiter_asteroid = %.5ei %.5ej km/s^2\n', dominant_jupiter_asteroid(1), dominant_jupiter_asteroid(2));
fprintf('mag_dominant_jupiter_asteroid = %.5e km/s^2\n', norm(dominant_jupiter_asteroid));

% perturbing due to sun
direct_sun = Gm_sun*r_asteroid_sun/(norm(r_asteroid_sun))^3;
indirect_sun = Gm_sun*r_jupiter_sun/(norm(r_jupiter_sun))^3;
net_pert_sun = direct_sun - indirect_sun;
fprintf('direct_sun = %.5ei %.5ej km/s^2\n', direct_sun(1), direct_sun(2));
fprintf('mag_direct_sun = %.5e km/s^2\n', norm(direct_sun));
fprintf('indirect_sun = %.5ei %.5ej km/s^2\n', indirect_sun(1), indirect_sun(2));
fprintf('mag_indirect_sun = %.5e km/s^2\n', norm(indirect_sun));
fprintf('net_pert_sun = %.5ei %.5ej km/s^2\n\n', net_pert_sun(1), net_pert_sun(2));
fprintf('mag_net_pert_sun = %.5e km/s^2\n\n', norm(net_pert_sun));

total_net_jupiter_asteroid = dominant_jupiter_asteroid + net_pert_sun;
fprintf('total net jupiter asteroid = %.5ei %.5ej km/s^2\n', total_net_jupiter_asteroid(1), total_net_jupiter_asteroid(2));
fprintf('mag total net jupiter asteroid = %.5e km/s^2\n\n', norm(total_net_jupiter_asteroid));