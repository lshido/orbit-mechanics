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

% Characteristics of hyperbola for jupiter flyby
specific_energy_flyby = v_inf_arr^2/2;
fprintf("specific_energy_flyby: %.4e km^2/s^2\n", specific_energy_flyby);

a_hyp_flyby = Gm_jupiter/(2*specific_energy_flyby);
fprintf("a_hyp_flyby: %.4e km\n", a_hyp_flyby);

e_hyp_flyby = (rp_jupiter/a_hyp_flyby) + 1;
fprintf("e_hyp_flyby: %.4e km\n", e_hyp_flyby);

% calculate flyby angle
turn_angle = asind(1/e_hyp_flyby);
fprintf("turn_angle: %.4e deg\n", turn_angle);
FBA = 2*asind(1/e_hyp_flyby);
fprintf("FBA: %.4e deg\n", FBA);

rp_and_a = rp_jupiter + a_hyp_flyby;
fprintf("rp_and_a: %.4e km\n", rp_and_a);

p_hyp_flyby = a_hyp_flyby*((e_hyp_flyby^2) - 1);
fprintf("p_hyp_flyby: %.4e km\n", p_hyp_flyby);

L = sqrt((rp_and_a^2) + (p_hyp_flyby^2));
fprintf("L: %.4e km\n", L);

kappa_angle = acosd(rp_and_a/L);
fprintf("kappa_angle: %.4e deg\n", kappa_angle);
kappa_angle_check = asind(p_hyp_flyby/L);
fprintf("kappa_angle_check: %.4e deg\n", kappa_angle_check);

FBA_check = 180 - (2*kappa_angle);
fprintf("FBA_check: %.4e deg\n", FBA_check);

b_hyp_flyby = a_hyp_flyby*sqrt((e_hyp_flyby^2)-1);
fprintf("b_hyp_flyby: %.4e deg\n", b_hyp_flyby);
kappa_angle_check_2 = asind(b_hyp_flyby/rp_and_a);
fprintf("kappa_angle_check_2: %.4e deg\n", kappa_angle_check_2);

FBA_check_2 = 180 - (2*kappa_angle_check_2);
fprintf("FBA_check_2: %.4e deg\n", FBA_check_2);

alpha_angle_flyby = (180-FBA)/2;
fprintf("alpha_angle_flyby: %.4e deg\n", alpha_angle_flyby);

delta_v_flyby = 2*v_inf_arr*cosd(alpha_angle_flyby);
fprintf("delta_v_flyby: %.4e km/s\n", delta_v_flyby);

v_plus_flyby = sqrt((va^2)+(delta_v_flyby^2)-2*va*delta_v_flyby*cosd(180-alpha_angle_flyby));
fprintf("v_plus_flyby: %.4e km/s\n", v_plus_flyby);

gamma_plus_flyby = asind((delta_v_flyby*sind(180-alpha_angle_flyby))/v_plus_flyby);
fprintf("gamma_plus_flyby: %.4e km/s\n", gamma_plus_flyby);

TA_plus_flyby = atand(((a_jupiter*v_plus_flyby^2/Gm_sun)*cosd(gamma_plus_flyby)*sind(gamma_plus_flyby))/((a_jupiter*v_plus_flyby^2/Gm_sun)*(cosd(gamma_plus_flyby))^2)+1);
fprintf("TA_plus_flyby: %.4e km/s\n", TA_plus_flyby);

specific_energy_flyby_new = ((v_plus_flyby^2)/2)-(Gm_sun/a_jupiter);
fprintf("specific_energy_flyby_new: %.4e km/s\n", specific_energy_flyby_new);

a_flyby_new = -Gm_sun/(2*specific_energy_flyby_new);
fprintf("a_flyby_new: %.4e km\n", a_flyby_new);

v_plus_flyby_rhat = v_plus_flyby*cosd(gamma_plus_flyby);
fprintf("v_plus_flyby_rhat: %.4e km/s\n", v_plus_flyby_rhat);

v_plus_flyby_thetahat = v_plus_flyby*sind(gamma_plus_flyby);
fprintf("v_plus_flyby_thetahat: %.4e km/s\n", v_plus_flyby_thetahat);

h = a_jupiter*v_plus_flyby_thetahat;
fprintf("h: %.4e km^2/s\n", h);

p_flyby_new = h^2/Gm_sun;
fprintf("p_flyby_new: %.4e km\n", p_flyby_new);

e_flyby_new = sqrt(1-(p_flyby_new/a_flyby_new));
fprintf("e_flyby_new: %.4e\n", e_flyby_new);

ra_flyby_new = a_flyby_new*(1+e_flyby_new);
rp_flyby_new = a_flyby_new*(1-e_flyby_new);
fprintf("ra_flyby_new: %.4e km\n", ra_flyby_new);
fprintf("rp_flyby_new: %.4e km\n", rp_flyby_new);

period_flyby_new = 2*pi*sqrt((a_flyby_new^3)/Gm_sun);
fprintf("period_flyby_new: %.4e sec\n", period_flyby_new);
fprintf("period_flyby_new: %.4e years\n", period_flyby_new/3600/24/365);

delta_small_omega = -TA_plus_flyby + 180;
fprintf("delta_small_omega: %.4e deg\n", delta_small_omega);

specific_energy_flyby_old = ((va^2)/2) - (Gm_sun/a_jupiter);
fprintf("specific_energy_flyby_old: %.4e km^2/s^2\n", specific_energy_flyby_old);