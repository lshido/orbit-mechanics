mu_mars = 42828.314258067;
R_mars = 3397;
a_phobos = 9376;
a_deimos = 23458;
TA = 240;

fprintf("-------Problem 1a----------\n")
r1 = a_phobos;
r2 = a_deimos;
phi_angle = 360-TA;
c = sqrt((r1^2) + (r2^2) - (2*r1*r2*cosd(phi_angle)));
s = (r1 + r2 + c)/2;
TOF_par_type2 = 1/3*sqrt(2/mu_mars)*((s^(3/2)) + ((s-c)^(3/2)));
a_min = s/2;
alpha_min = 2*asind(sqrt(s/(2*a_min)));
beta_min = 2*asind(sqrt((s-c)/(2*a_min)));
TOF_min = sqrt(a_min^3/mu_mars)*((deg2rad(alpha_min) - sind(alpha_min))+(deg2rad(beta_min) - sind(beta_min)));
fprintf("c: %.4e km\n", c)
fprintf("s: %.4e km\n", s)
fprintf("TOF_par_type2: %.4e sec\n", TOF_par_type2)
fprintf("TOF_par_type2: %.4e hrs\n", TOF_par_type2/3600)
fprintf("s: %.4e km\n", s)
fprintf("a_min: %.4e km\n", a_min)
fprintf("alpha_min: %.4e rad\n", deg2rad(alpha_min))
fprintf("alpha_min: %.4e deg\n", alpha_min)
fprintf("beta_min: %.4e rad\n", deg2rad(beta_min))
fprintf("beta_min: %.4e deg\n", beta_min)
fprintf("TOF_min: %.4e sec\n", TOF_min)
fprintf("TOF_min: %.4e hrs\n", TOF_min/3600)

fprintf("-------Problem 1a: Iterate on 'a'----------\n")
a_guess = a_min;
step_size = .01;
TOF_desired = 15*3600;
accuracy = 1;
TOF_error = 2;

while abs(TOF_error) > accuracy

    if TOF_error > 0
        a_guess = a_guess - step_size;
    else
        a_guess = a_guess + step_size;
    end

    alpha_guess = 2*asind(sqrt(s/(2*a_guess)));
    beta_guess = 2*asind(sqrt((s-c)/(2*a_guess)));
    TOF_guess = sqrt(a_guess^3/mu_mars)*((2*pi) - (deg2rad(alpha_guess) - sind(alpha_guess)) + (deg2rad(beta_guess) - sind(beta_guess)));
    TOF_error = TOF_guess - TOF_desired;

end
fprintf("alpha_guess: %.4e rad\n", deg2rad(alpha_guess));
fprintf("alpha_guess: %.4e deg\n", alpha_guess);
fprintf("beta_guess: %.4e rad\n", deg2rad(beta_guess));
fprintf("beta_guess: %.4e deg\n,", beta_guess);
fprintf("a_guess: %.4e km,", a_guess);
fprintf("TOF_error: %.4e\n", TOF_error);

fprintf("-------Problem 1a: Calculate 'p's----------\n")
p_plus = ((4*a_guess*(s-r1)*(s-r2))/(c^2)) * (sind((alpha_guess + beta_guess)/2))^2;
p_minus = ((4*a_guess*(s-r1)*(s-r2))/(c^2)) * (sind((alpha_guess - beta_guess)/2))^2;
p = p_plus;
e = sqrt(1- (p/a_guess));
energy = -mu_mars/(2*a_guess);
rp = a_guess*(1-e);
ra = a_guess*(1+e);
fprintf("energy: %.4e km^2/s^2\n", energy);
fprintf("p_plus: %.4e km\n", p_plus);
fprintf("p_minus: %.4e km\n", p_minus);
fprintf("p: %.4e km\n", p);
fprintf("e: %.4e\n", e);
fprintf("rp: %.4e km\n", rp);
fprintf("rp: %.4e R_mars\n", rp/R_mars);
fprintf("ra: %.4e km\n", ra);

fprintf("-------Problem 1a: Calculate velocities-----------\n")
v_dep = sqrt(2*(energy + (mu_mars/r1)));
v_arr = sqrt(2*(energy + (mu_mars/r2)));
fprintf("v_dep: %.4e km/s\n", v_dep);
fprintf("v_arr: %.4e km/s\n", v_arr);

fprintf("-------Problem 1a: Calculate true anomalies-----------\n")
theta_star_dep_1 = acosd((p-r1)/(r1*e));
theta_star_dep_2 = 360 - theta_star_dep_1;
theta_star_arr_1 = acosd((p-r2)/(r2*e));
theta_star_arr_2 = 360 - theta_star_arr_1;
fprintf("theta_star_dep_1: %.4e deg\n", theta_star_dep_1);
fprintf("theta_star_dep_2: %.4e deg\n", theta_star_dep_2);
fprintf("theta_star_arr_1: %.4e deg\n", theta_star_arr_1);
fprintf("theta_star_arr_2: %.4e deg\n", theta_star_arr_2);

fprintf("-------Problem 1a: Calculate FPAs-----------\n")
FPA_dep = acosd(sqrt(p*mu_mars)/(r1*v_dep));
FPA_arr = acosd(sqrt(p*mu_mars)/(r2*v_arr));
fprintf("FPA_dep: %.4e deg\n", FPA_dep);
fprintf("FPA_arr: %.4e deg\n", FPA_arr);




