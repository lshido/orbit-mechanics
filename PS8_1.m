R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
a_jupiter = 778279959; % km
a_Earth = 149597898; % km
R_jupiter = 71492; % km
Gm_jupiter = 126712767.8578; %km^3/sec^2
Gm_sun = 132712440017.99;


% Determine conditions at Earth before maneuver
% circular orbit
rp_earth = R_earth + 250;
vc_earth = sqrt(Gm_earth/rp_earth);
fprintf("rp_earth: %.4e km\n", rp_earth);
fprintf("vc_earth: %.4e km/s\n", vc_earth);

r_1 = a_Earth;
r_2 = a_jupiter;

% Determine transfer ellipse
r_1 = a_Earth;
r_2 = a_jupiter;
a_transfer = 1/2*(r_1 + r_2);
fprintf("a_transfer: %.4e km\n", a_transfer);
e_transfer = 1 - (r_1/a_transfer);
fprintf("e_transfer: %.4e \n", e_transfer);

% conditions after maneuver at 1
v_earth = sqrt(Gm_sun/a_Earth);
fprintf("v_earth: %.4e km/sec\n", v_earth);
v_inf_dep = sqrt(2*((Gm_earth/r_1)-(Gm_earth/(2*a_transfer))));
fprintf("v_inf_dep: %.4e km/sec\n", v_inf_dep);
v_dep = v_earth + v_inf_dep;
fprintf("v_dep: %.4e km/sec\n", v_dep);

% delta v dep
delta_v_dep = sqrt(v_inf_dep^2 + (2*Gm_earth/rp_earth))-vc_earth;
fprintf("delta_v_dep: %.4e km/sec\n", delta_v_dep);

% conditions before maneuver at 2
v_jupiter = sqrt(Gm_sun/a_jupiter);
fprintf("v_jupiter: %.4e km/sec\n", v_jupiter);
v_inf_arr = sqrt(2*((Gm_earth/r_2)-(Gm_earth/(2*a_transfer))));
fprintf("v_inf_arr: %.4e km/sec\n", v_inf_arr);
v_arr = v_jupiter - v_inf_arr;
fprintf("v_arr: %.4e km/sec\n", v_arr);

% conditions required after maneuver in final orbit

rp_jupiter = 2.8*R_jupiter;
vc_jupiter = sqrt(Gm_earth/rp_jupiter);
fprintf("rp_jupiter: %.4e km\n", rp_jupiter);
fprintf("vc_jupiter: %.4e km/s\n", vc_jupiter);

% delta v arrival
delta_v_arr = sqrt(v_inf_dep^2 + (2*Gm_jupiter/rp_jupiter))-vc_jupiter;
fprintf("delta_v_arr: %.4e km/sec\n", delta_v_arr);

