Gm_sun = 132712440017.99;
R_venus = 6051.9;

earth_launch = 2453221.0;
first_venus_flyby = 2454033.0;
first_mercury_flyby = 2454480.0;
mercury_orbit_insertion = 2455639.0;

fprintf("-----------Problem 2a: TOFs--------\n")
earth_launch_to_first_venus_flyby = first_venus_flyby - earth_launch;
earth_launch_to_first_mercury_flyby = first_mercury_flyby - earth_launch;
earth_launch_to_mercury_orbit_insertion = mercury_orbit_insertion - earth_launch;
fprintf("earth_launch_to_first_venus_flyby: %.4e days\n", earth_launch_to_first_venus_flyby);
fprintf("earth_launch_to_first_venus_flyby: %.4e years\n", earth_launch_to_first_venus_flyby/365);
fprintf("earth_launch_to_first_mercury_flyby: %.4e days\n", earth_launch_to_first_mercury_flyby);
fprintf("earth_launch_to_first_mercury_flyby: %.4e years\n", earth_launch_to_first_mercury_flyby/365);
fprintf("earth_launch_to_mercury_orbit_insertion: %.4e days\n", earth_launch_to_mercury_orbit_insertion);
fprintf("earth_launch_to_mercury_orbit_insertion: %.4e years\n", earth_launch_to_mercury_orbit_insertion/365);

fprintf("-----------Problem 2b: Compute TA at venus encounter--------\n")
rp_Hohmann = 7.4800e7;
ra_Hohmann = 1.4960e8;
r_encounter = 1.0821e8;
a_Hohmann = (rp_Hohmann + ra_Hohmann)/2;
e_Hohmann = 1 - (rp_Hohmann/a_Hohmann);
p_Hohmann = a_Hohmann*(1-(e_Hohmann^2));
TA_encounter = acosd((p_Hohmann/(r_encounter*e_Hohmann))-(1/e_Hohmann));
EA_encounter = pi + 2*atan(sqrt((1-e_Hohmann)/(1+e_Hohmann))*tand(TA_encounter/2));
time_to_venus = sqrt((a_Hohmann^3)/Gm_sun)*(EA_encounter - (e_Hohmann*sin(EA_encounter)));
time_to_earth = sqrt((a_Hohmann^3)/Gm_sun)*pi;
TOF_earth_to_venus = time_to_venus - time_to_earth;
fprintf("rp_Hohmann: %.4e km\n", rp_Hohmann);
fprintf("ra_Hohmann: %.4e km\n", ra_Hohmann);
fprintf("a_Hohmann: %.4e km\n", a_Hohmann);
fprintf("e_Hohmann: %.4e\n", e_Hohmann);
fprintf("p_Hohmann: %.4e\n", p_Hohmann);
fprintf("TA_encounter: %.4e deg\n", TA_encounter);
fprintf("EA_encounter: %.4e rad\n", EA_encounter);
fprintf("EA_encounter: %.4e deg\n", rad2deg(EA_encounter));
fprintf("time_to_venus: %.4e sec\n", time_to_venus);
fprintf("time_to_earth: %.4e sec\n", time_to_earth);
fprintf("TOF_earth_to_venus: %.4e sec\n", TOF_earth_to_venus);
fprintf("TOF_earth_to_venus: %.4e days\n", TOF_earth_to_venus/3600/24);

fprintf("-----------Problem 2c: Computer post-Venus Encounter--------\n")
rp_venus = R_venus + 2990;
fprintf("rp_venus: %.4e km\n", rp_venus);
fprintf("rp_venus radii: %.4e km\n", rp_venus/R_venus);
