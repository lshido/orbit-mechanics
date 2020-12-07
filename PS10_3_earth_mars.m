mu_sun = 132712440017.99;
mu_earth = 398600.4415;
mu_mars = 42828.314258067;

R_sun = 695990;
R_earth = 6378.1363;
R_mars = 3397;

a_earth = 149597898;
a_mars = 227944135;
TA = 120;

fprintf("-------Problem 3a----------\n")
r1 = a_earth;
r2 = a_mars;
phi_angle = TA;
c = sqrt((r1^2) + (r2^2) - (2*r1*r2*cosd(phi_angle)));
s = (r1 + r2 + c)/2;
TOF_par_type1 = 1/3*sqrt(2/mu_sun)*((s^(3/2)) - ((s-c)^(3/2)));
fprintf("c: %.4e km\n", c)
fprintf("s: %.4e km\n", s)
fprintf("TOF_par_type1: %.4e sec\n", TOF_par_type1)
fprintf("TOF_par_type1: %.4e hrs\n", TOF_par_type1/3600)
fprintf("TOF_par_type1: %.4e days\n", TOF_par_type1/3600/24)

fprintf("-------Problem 3b: Iterate on 'a'----------\n")
a_guess = 0.1;
step_size = 1000;
TOF_desired = 92*24*3600;
accuracy = 1;
TOF_error = 11;

while 1

    if abs(TOF_error) <= accuracy
        break
    elseif TOF_error > 0
        a_guess = a_guess + step_size;
    elseif TOF_error < 0
        a_guess = a_guess - step_size;
    end

    alpha0_guess = 2*asinh(sqrt(s/(2*abs(a_guess))));
    beta0_guess = 2*asinh(sqrt((s-c)/(2*abs(a_guess))));
    TOF_guess = sqrt(abs(a_guess)^3/mu_sun)*((sinh(alpha0_guess) - alpha0_guess) - (sinh(beta0_guess) - beta0_guess));
    TOF_error = TOF_guess - TOF_desired;
end
alpha_guess = alpha0_guess;
beta_guess = beta0_guess;
a = abs(a_guess);
fprintf("a_guess: %.4e km,", a_guess);
fprintf("TOF_error: %.4e\n", TOF_error);
fprintf("a: %.4e km\n", a);
fprintf("alpha0_guess: %.4e rad\n", alpha0_guess);
fprintf("alpha0_guess: %.4e deg\n", rad2deg(alpha0_guess));
fprintf("beta0_guess: %.4e rad\n", beta0_guess);
fprintf("beta0_guess: %.4e deg\n", rad2deg(beta0_guess));
fprintf("alpha_guess: %.4e deg\n", alpha_guess);
fprintf("beta_guess: %.4e deg\n", beta_guess);

fprintf("-------Problem 3b: Calculate 'p's----------\n")
p_plus = ((4*a*(s-r1)*(s-r2))/(c^2)) * (sinh((alpha_guess + beta_guess)/2))^2;
p_minus = ((4*a*(s-r1)*(s-r2))/(c^2)) * (sinh((alpha_guess - beta_guess)/2))^2;
p = p_plus;
e = sqrt(1 + (p/a));
energy = mu_sun/(2*a);
rp = a*(e-1);
fprintf("energy: %.4e km^2/s^2\n", energy);
fprintf("p_plus: %.4e km\n", p_plus);
fprintf("p_minus: %.4e km\n", p_minus);
fprintf("p: %.4e km\n", p);
fprintf("e: %.4e\n", e);
fprintf("rp: %.4e km\n", rp);
fprintf("rp: %.4e R_sun\n", rp/R_sun);

fprintf("-------Problem 3b: Calculate velocities-----------\n")
v_dep = sqrt(2*(energy + (mu_sun/r1)));
v_arr = sqrt(2*(energy + (mu_sun/r2)));
fprintf("v_dep: %.4e km/s\n", v_dep);
fprintf("v_arr: %.4e km/s\n", v_arr);

fprintf("-------Problem 3b: Calculate true anomalies-----------\n")
theta_star_dep_1 = acosd((p-r1)/(r1*e));
theta_star_dep_2 = 360 - theta_star_dep_1;
theta_star_arr_1 = acosd((p-r2)/(r2*e));
theta_star_arr_2 = 360 - theta_star_arr_1;
fprintf("theta_star_dep_1: %.4e deg\n", theta_star_dep_1);
fprintf("theta_star_dep_2: %.4e deg\n", theta_star_dep_2);
fprintf("theta_star_arr_1: %.4e deg\n", theta_star_arr_1);
fprintf("theta_star_arr_2: %.4e deg\n", theta_star_arr_2);

fprintf("-------Problem 23b: Calculate FPAs-----------\n")
FPA_dep = acosd(sqrt(p*mu_sun)/(r1*v_dep));
FPA_arr = acosd(sqrt(p*mu_sun)/(r2*v_arr));
fprintf("FPA_dep: %.4e deg\n", FPA_dep);
fprintf("FPA_arr: %.4e deg\n", FPA_arr);