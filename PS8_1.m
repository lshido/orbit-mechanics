R_earth = 6378.1363; % [km]
Gm_earth = 398600.4415;
a_jupiter = 778279959; % km
a_Earth = 149597898; % km
R_jupiter = 71492; % km
Gm_jupiter = 126712767.8578; %km^3/sec^2
Gm_sun = 132712440017.99;


rp_earth = R_earth + 250;
rp_jupiter = 2.8*R_jupiter;

% circular orbit of Earth
v_earth = sqrt(Gm_sun/a_Earth);
fprintf("v_earth: %.4e km/s\n", v_earth);
v_earth_vector = [ v_earth 0 ];

% Determine transfer ellipse
r_1 = a_Earth;
r_2 = a_jupiter;
a_transfer = 1/2*(r_1 + r_2);
fprintf("a_transfer: %.4e km\n", a_transfer);
e_transfer = 1 - (r_1/a_transfer);
fprintf("e_transfer: %.4e \n", e_transfer);

% Calculate velocity at ellipse periapsis
vp = sqrt((2*Gm_sun/r_1)-(Gm_sun/a_transfer));
fprintf("vp: %.4e km/s\n", vp);
vp_vector = [ vp 0 ];

% calculate v_infinity leaving earth
v_inf_dep_vector = vp_vector - v_earth_vector;
v_inf_dep = norm(v_inf_dep_vector);
fprintf("mag of v_inf_dep_vector: %.4e km/s\n", v_inf_dep);

% calculate velocity in parking orbit
v_parking_earth = sqrt(Gm_earth/rp_earth);
fprintf("v_parking_earth: %.4e km/s\n", v_parking_earth);

% calculate delta_v_dep
delta_v_dep = sqrt((v_inf_dep^2) + (2*Gm_earth/rp_earth)) - v_parking_earth;
fprintf("delta_v_dep: %.4e km/s\n", delta_v_dep);

% circular orbit of Jupiter
v_jupiter = sqrt(Gm_sun/a_jupiter);
fprintf("v_jupiter: %.4e km/s\n", v_jupiter);
v_jupiter_vector = [ v_jupiter 0 ];

% Calculate velocity at ellipse apoapsis
va = sqrt((2*Gm_sun/r_2)-(Gm_sun/a_transfer));
fprintf("va: %.4e km/s\n", va);
va_vector = [ va 0 ];

% calculate v_infinity leaving earth
v_inf_arr_vector = va_vector - v_jupiter_vector;
v_inf_arr = norm(v_inf_arr_vector);
fprintf("v_inf_arr_vector: %.4e km/s %.4e km/s\n", v_inf_arr_vector);
fprintf("mag of v_inf_arr_vector: %.4e km/s\n", v_inf_arr);

% calculate velocity in parking orbit
v_parking_jupiter = sqrt(Gm_jupiter/rp_jupiter);
fprintf("v_parking_jupiter: %.4e km/s\n", v_parking_jupiter);

% calculate delta_v_dep
delta_v_arr = sqrt((v_inf_arr^2) + (2*Gm_jupiter/rp_jupiter)) - v_parking_jupiter;
fprintf("delta_v_arr: %.4e km/s\n", delta_v_arr);

% calculate delta_v_total
delta_v_total = delta_v_dep + delta_v_arr;
fprintf("delta_v_total: %.4e km/s\n", delta_v_total);

% calculate TOF
TOF_hohmann = pi*sqrt(a_transfer^3/Gm_sun);
fprintf("TOF_hohmann in sec: %.4e sec\n", TOF_hohmann);
fprintf("TOF_hohmann in days: %.4e days\n", TOF_hohmann/3600/24);
fprintf("TOF_hohmann in years: %.4e years\n", TOF_hohmann/3600/24/365);

% a_jupiter
e_jupiter = 0.90;
a_eccentric_jupiter = rp_jupiter/(1-e_jupiter);
fprintf("a_eccentric_jupiter: %.4e km\n", a_eccentric_jupiter);

% Calculate velocity at perijove
vp_eccentric_jupiter = sqrt((2*Gm_jupiter/rp_jupiter)-(Gm_jupiter/a_eccentric_jupiter));
fprintf("vp_eccentric_jupiter: %.4e km/s\n", vp_eccentric_jupiter);

% calculate delta_v_arr_eccentric
delta_v_arr_eccentric = sqrt((v_inf_arr^2) + (2*Gm_jupiter/rp_jupiter)) - vp_eccentric_jupiter;
fprintf("delta_v_arr_eccentric: %.4e km/s\n", delta_v_arr_eccentric);

% delta v arr change
delta_v_arr_change = delta_v_arr_eccentric - delta_v_arr;
fprintf("delta_v_arr_change: %.4e km/s\n", delta_v_arr_change);

% delta v total new = 
delta_v_total_new = delta_v_total + delta_v_arr_change;
fprintf("delta_v_total_new: %.4e km/s\n", delta_v_total_new);
