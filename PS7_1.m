R_mars = 3397; % km
Gm_mars = 42828.314258067; % km^3/sec^2

a_old = 6*R_mars;
fprintf("Semi major axis before maneuver: %.4e km", a_old);

specific_energy_old = -Gm_mars/(2*a_old);
fprintf("specific energy before maneuver: %.4e km^2/s^2", specific_energy_old);
