Gm_sun = 132712440017.99;
Gm_earth = 398600.4415;
Gm_jupiter = 126712767.8578;
R_earth = 6378.1363; % [km]

a_Earth = 149597898; % km
%a_Earth = 2*R_earth;
a_jupiter = 778279959; % km
%a_jupiter = 4*R_earth;
e_jupiter = 0.04853590;

launch = 2455779.0;
arrival = 2457575.0;
EFB = 2456575.0;

TOF_days = arrival - launch;
fprintf("TOF in days: %.8e\n", TOF_days);
fprintf("TOF in years: %.8e\n", TOF_days/365);
TOF_adjusted = arrival - EFB;
fprintf("TOF_adjusted in days: %.8e\n", TOF_adjusted);
fprintf("TOF_adjusted in years: %.8e\n", TOF_adjusted/365);
r_1 = a_Earth;
r_2 = a_jupiter;

% establish current orbit
e = 0;

% conditions before maneuver at 1
r_1_before = a_Earth;
FPA_1_before = 0; % deg
%mu_1 = Gm_sun + Gm_earth;
mu_1 = Gm_sun;
specific_energy_1_before = -mu_1/(2*a_Earth);
fprintf("specific_energy_1_before: %.4e\n", specific_energy_1_before);

% calc v_1_before
v_1_before = sqrt(2*(specific_energy_1_before+(mu_1/r_1_before)));
fprintf("v_1_before: %.4e km/sec\n", v_1_before);

% transfer ellipse
a_transfer = 1/2*(r_1 + r_2);
fprintf("a_transfer: %.4e km\n", a_transfer);
e_transfer = 1 - (r_1/a_transfer);
fprintf("e_transfer: %.4e \n", e_transfer);

% conditions after maneuver at 1
v_1_after = sqrt(2*((mu_1/r_1)-(mu_1/(2*a_transfer))));
fprintf("v_1_after: %.4e km/sec\n", v_1_after);

% calculate delta v 1
delta_v_1 = v_1_after - v_1_before;
fprintf("delta_v_1: %.4e km/sec\n", delta_v_1);

% conditions before maneuver at 2
v_2_before = sqrt(2*((mu_1/r_2)-(mu_1/(2*a_transfer))));
fprintf("v_2_before: %.4e km/sec\n", v_2_before);

% conditions required after maneuver in final orbit
v_2_after = sqrt(mu_1/r_2);
fprintf("v_2_after: %.4e km/sec\n", v_2_after);
delta_v_2 = v_2_after - v_2_before;
fprintf("delta_v_2: %.4e km/sec\n", delta_v_2);

% calculate total delta_v
delta_v_total = delta_v_1 + delta_v_2;
fprintf("delta_v_total: %.4e km/sec\n", delta_v_total);

% calculate TOF
TOF_hohmann = pi*sqrt(a_transfer^3/mu_1);
fprintf("TOF_hohmann in sec: %.8e sec\n", TOF_hohmann);
fprintf("TOF_hohmann in days: %.8e days\n", TOF_hohmann/3600/24);
fprintf("TOF_hohmann in years: %.8e years\n", TOF_hohmann/3600/24/365);

% calculate phase angle phi
n_2 = sqrt(mu_1/a_jupiter^3);
fprintf("n_2: %.4e rad/sec\n", n_2);
phi = 180 - rad2deg(n_2*TOF_hohmann);
fprintf("phi: %.4e deg\n", phi);

n_1 = sqrt(mu_1/a_Earth^3);
fprintf("n_1: %.4e rad/sec\n", n_1);
synodic = 2*pi/(n_1 - n_2);
fprintf("synodic in sec: %.4e sec\n", synodic);
fprintf("synodic in days: %.4e days\n", synodic/3600/24);
fprintf("synodic in years: %.4e years\n", synodic/3600/24/365);

% assume eccentric jupiter orbit
r_p_jupiter = a_jupiter*(1-e_jupiter);
fprintf("r_p_jupiter: %.4e km\n", r_p_jupiter);

% determine transfer ellipse
% transfer ellipse
a_transfer_eccentric = 1/2*(r_1 + r_p_jupiter);
fprintf("a_transfer_eccentric: %.4e km\n", a_transfer_eccentric);
e_transfer_eccentric = 1 - (r_1/a_transfer_eccentric);
fprintf("e_transfer_eccentric: %.4e \n", e_transfer_eccentric);

% conditions after maneuver at 1
v_1_after_eccentric = sqrt(2*((mu_1/r_1)-(mu_1/(2*a_transfer_eccentric))));
fprintf("v_1_after_eccentric: %.4e km/sec\n", v_1_after_eccentric);

% calculate delta v 1
delta_v_1_eccentric = v_1_after_eccentric - v_1_before;
fprintf("delta_v_1_eccentric: %.4e km/sec\n", delta_v_1_eccentric);

% conditions before maneuver at 2
v_2_before_eccentric = sqrt(2*((mu_1/r_p_jupiter)-(mu_1/(2*a_transfer_eccentric))));
fprintf("v_2_before_eccentric: %.4e km/sec\n", v_2_before_eccentric);

% conditions required after maneuver in final orbit
v_2_after_eccentric = sqrt(2*((mu_1/r_p_jupiter)-(mu_1/(2*a_jupiter))));
fprintf("v_2_after_eccentric: %.4e km/sec\n", v_2_after_eccentric);

delta_v_2_eccentric = v_2_after_eccentric - v_2_before_eccentric;
fprintf("delta_v_2_eccentric: %.4e km/sec\n", delta_v_2_eccentric);

% calculate total delta_v
delta_v_total_eccentric = delta_v_1_eccentric + delta_v_2_eccentric;
fprintf("delta_v_total_eccentric: %.4e km/sec\n", delta_v_total_eccentric);

%calc TOF hohmann with eccentric orbit
TOF_hohmann_eccentric = pi*sqrt(a_transfer_eccentric^3/mu_1);
fprintf("TOF_hohmann_eccentric in sec: %.8e sec\n", TOF_hohmann_eccentric);
fprintf("TOF_hohmann_eccentric in days: %.8e days\n", TOF_hohmann_eccentric/3600/24);
fprintf("TOF_hohmann_eccentric in years: %.8e years\n", TOF_hohmann_eccentric/3600/24/365);

