mu_mars = 42828.37;
a_phobos = 9376;
a_deimos = 23458;
TA = 240;

fprintf("-------Problem 1a----------\n")
r1 = a_phobos;
r2 = a_deimos;
phi_angle = 360-TA;
c = sqrt(r1^2 + r2^2 - 2*r1*r2*cosd(phi_angle));
s = (r1 + r2 + c)/2;
TOF_par_type2 = 1/3*sqrt(2/mu_mars)*(s^(3/2) + (s-c)^(3/2));
a_min = s/2;
alpha_min = 2*asin(s/(2*a_min));
beta_min = 2*asin((s-c)/(2*a_min));
TOF_min = sqrt(a_min^3/mu_mars)*((alpha_min - sin(alpha_min))-(beta_min - sin(beta_min)));
fprintf("c: %.4e km\n", c)
fprintf("s: %.4e km\n", s)
fprintf("TOF_par_type2: %.4e sec\n", TOF_par_type2)
fprintf("s: %.4e km\n", s)
fprintf("a_min: %.4e km\n", a_min)
fprintf("alpha_min: %.4e rad\n", alpha_min)
fprintf("alpha_min: %.4e deg\n", rad2deg(alpha_min))
fprintf("beta_min: %.4e rad\n", beta_min)
fprintf("beta_min: %.4e deg\n", rad2deg(beta_min))
fprintf("TOF_min: %.4e sec\n", TOF_min)
fprintf("TOF_min: %.4e hrs\n", TOF_min/3600)

fprintf("-------Problem 1a: Iterate on 'a'----------\n")
a_guess = a_min;
step_size = .1;
TOF_desired = 15*3600;
accuracy = 1;
TOF_error = 100;

while abs(TOF_error) > accuracy

    if TOF_error > 0
        a_guess = a_guess - step_size;
    else
        a_guess = a_guess + step_size;
    end

    alpha_guess = 2*asin(s/(2*a_guess));
    beta_guess = 2*asin((s-c)/(2*a_guess));
    TOF_guess = sqrt(a_guess^3/mu_mars)*(2*pi - (alpha_guess - sin(alpha_guess)) - (beta_guess - sin(beta_guess)));
    TOF_error = TOF_guess - TOF_desired;

end
fprintf("alpha_guess: %.4e rad\n", alpha_guess);
fprintf("alpha_guess: %.4e deg\n", rad2deg(alpha_guess));
fprintf("beta_guess: %.4e rad\n", beta_guess);
fprintf("beta_guess: %.4e deg\n,", rad2deg(beta_guess));
fprintf("a_guess: %.4e km,", a_guess);
fprintf("TOF_error: %.4e\n", TOF_error);

fprintf("-------Problem 1a: Calculate 'p's----------\n")
p_plus = (4*a_guess*(s-r1)*(s-r2))/c^2 * (sin((alpha_guess + beta_guess)/2))^2;
p_minus = (4*a_guess*(s-r1)*(s-r2))/c^2 * (sin((alpha_guess - beta_guess)/2))^2;
p = p_plus;
e = sqrt(1- (p/a_guess));
fprintf("p_plus: %.4e km\n", p_plus);
fprintf("p_minus: %.4e km\n", p_minus);
fprintf("p: %.4e km\n", p);
fprintf("e: %.4e\n", e);




