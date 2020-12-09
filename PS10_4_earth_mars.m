mu_sun = 132712440017.99;
mu_earth = 398600.4415;
mu_mars = 42828.314258067;

R_sun = 695990;
R_earth = 6378.1363;
R_mars = 3397;

a_earth = 149597898;
a_mars = 227944135;
TA = 60;

fprintf("-------Problem 4a----------\n")
r1 = a_earth;
r2 = a_mars;
phi_angle = TA;
c = sqrt((r1^2) + (r2^2) - (2*r1*r2*cosd(phi_angle)));
s = (r1 + r2 + c)/2;
TOF_par_type1 = 1/3*sqrt(2/mu_sun)*((s^(3/2)) - ((s-c)^(3/2)));
a_min = s/2;
energy_min = -mu_sun/(2*a_min);
alpha_min = 2*asind(sqrt(s/(2*a_min)));
beta_min = 2*asind(sqrt((s-c)/(2*a_min)));
TOF_min = sqrt(a_min^3/mu_sun)*((deg2rad(alpha_min) - sind(alpha_min))-(deg2rad(beta_min) - sind(beta_min)));
fprintf("c: %.4e km\n", c)
fprintf("s: %.4e km\n", s)
fprintf("TOF_par_type1: %.4e sec\n", TOF_par_type1)
fprintf("TOF_par_type1: %.4e hrs\n", TOF_par_type1/3600)
fprintf("TOF_par_type1: %.4e days\n", TOF_par_type1/3600/24)
fprintf("a_min: %.4e km\n", a_min)
fprintf("alpha_min: %.4e rad\n", deg2rad(alpha_min))
fprintf("alpha_min: %.4e deg\n", alpha_min)
fprintf("beta_min: %.4e rad\n", deg2rad(beta_min))
fprintf("beta_min: %.4e deg\n", beta_min)
fprintf("TOF_min: %.4e sec\n", TOF_min)
fprintf("TOF_min: %.4e hrs\n", TOF_min/3600)
fprintf("TOF_min: %.4e days\n", TOF_min/3600/24)

fprintf("-------Problem 4b: Iterate on 'a'----------\n")
a_guess = a_min;
step_size = 10;
TOF_desired = 92*24*3600;
accuracy = 1;
TOF_error = 11;

while 1

    if TOF_error <= accuracy
        break
    elseif TOF_error > 0
        a_guess = a_guess + step_size;
    elseif TOF_error < 0
        a_guess = a_guess - step_size;
    end

    alpha0_guess = 2*asind(sqrt(s/(2*a_guess)));
    beta0_guess = 2*asind(sqrt((s-c)/(2*a_guess)));
    TOF_guess = sqrt(a_guess^3/mu_sun)*((deg2rad(alpha0_guess) - sind(alpha0_guess)) - (deg2rad(beta0_guess) - sind(beta0_guess)));
    TOF_error = TOF_guess - TOF_desired;
end
alpha_guess = alpha0_guess;
beta_guess = beta0_guess;
fprintf("a_guess: %.4e km,", a_guess);
fprintf("TOF_error: %.4e\n", TOF_error);
fprintf("alpha0_guess: %.4e rad\n", deg2rad(alpha0_guess));
fprintf("alpha0_guess: %.4e deg\n", alpha0_guess);
fprintf("beta0_guess: %.4e rad\n", deg2rad(beta0_guess));
fprintf("beta0_guess: %.4e deg\n", beta0_guess);
fprintf("alpha_guess: %.4e deg\n", alpha_guess);
fprintf("beta_guess: %.4e deg\n", beta_guess);

fprintf("-------Problem 3b: Calculate 'p's----------\n")
p_plus = ((4*a_guess*(s-r1)*(s-r2))/(c^2)) * (sind((alpha_guess + beta_guess)/2))^2;
p_minus = ((4*a_guess*(s-r1)*(s-r2))/(c^2)) * (sind((alpha_guess - beta_guess)/2))^2;
p = p_plus;
e = sqrt(1- (p/a_guess));
energy = -mu_sun/(2*a_guess);
rp = a_guess*(1-e);
ra = a_guess*(1+e);
fprintf("energy: %.4e km^2/s^2\n", energy);
fprintf("p_plus: %.4e km\n", p_plus);
fprintf("p_minus: %.4e km\n", p_minus);
fprintf("p: %.4e km\n", p);
fprintf("e: %.4e\n", e);
fprintf("rp: %.4e km\n", rp);
fprintf("rp: %.4e R_sun\n", rp/R_sun);
fprintf("ra: %.4e km\n", ra);

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

fprintf("-------Problem 3b: Calculate FPAs-----------\n")
FPA_dep = acosd(sqrt(p*mu_sun)/(r1*v_dep));
FPA_arr = acosd(sqrt(p*mu_sun)/(r2*v_arr));
fprintf("FPA_dep: %.4e deg\n", FPA_dep);
fprintf("FPA_arr: %.4e deg\n", FPA_arr);

fprintf("-------Problem 3c: Calculate delta_v-----------\n")
v1 = sqrt(mu_sun/r1);
v2 = sqrt(mu_sun/r2);
delta_v1 = sqrt((v1^2)+(v_dep^2)-(2*v1*v_dep*cosd(FPA_dep)));
delta_v2 = sqrt((v2^2)+(v_arr^2)-(2*v2*v_arr*cosd(FPA_arr)));

kappa_angle_1 = asind((v_dep*sind(FPA_dep))/delta_v1);
kappa_angle_1_check = acosd(((v_dep^2) - (delta_v1^2) - (v1^2))/(-2*v1*delta_v1));
alpha_angle_1 = 180 - kappa_angle_1_check;
delta_v1_vector = [ delta_v1*cosd(alpha_angle_1) delta_v1*sind(alpha_angle_1) ];
kappa_angle_2 = asind((v2*sind(FPA_arr))/delta_v2);
kappa_angle_2_check = acosd(((v2^2) - (v_arr^2) - (delta_v2^2))/(-2*v_arr*delta_v2));
alpha_angle_2 = -(180 - kappa_angle_2_check);
delta_v2_vector = [ delta_v2*cosd(alpha_angle_2) delta_v2*sind(alpha_angle_2) ];

fprintf("v1: %.4e km/s\n", v1);
fprintf("v2: %.4e km/s\n", v2);
fprintf("delta_v1: %.4e km/s\n", delta_v1);
fprintf("delta_v2: %.4e km/s\n", delta_v2);

fprintf("kappa_angle_1: %.4e deg\n", kappa_angle_1);
fprintf("kappa_angle_1_check: %.4e deg\n", kappa_angle_1_check);
fprintf("alpha_angle_1: %.4e deg\n", alpha_angle_1);
fprintf("delta_v1_vector: %.4e V^ %.4e B^ km/s\n", delta_v1_vector);
fprintf("kappa_angle_2: %.4e deg\n", kappa_angle_2);
fprintf("kappa_angle_2_check: %.4e deg\n", kappa_angle_2_check);
fprintf("alpha_angle_2: %.4e deg\n", alpha_angle_2);
fprintf("delta_v2_vector: %.4e V^ %.4e B^ km/s\n", delta_v2_vector);

fprintf("-------Problem 2f: Consider Earth Local Field-----------\n")
r_pass_earth = 1000 + R_earth;
v_earth_vector = [ v1 0 ];
v_dep_vector = [ v_dep*cosd(FPA_dep) v_dep*sind(FPA_dep) ];
v_inf_earth_vector = v_dep_vector - v_earth_vector;
v_inf_earth = norm(v_inf_earth_vector);
energy_hyp_earth = v_inf_earth^2/2;
v_pass_earth = sqrt(2*(energy_hyp_earth + (mu_earth/r_pass_earth)));
fprintf("energy_hyp_earth: %.4e km/s\n", energy_hyp_earth);
fprintf("v_earth_vector: %.4e V^ %.4e B^ km/s\n", v_earth_vector);
fprintf("v_dep_vector: %.4e V^ %.4e B^ km/s\n", v_dep_vector);
fprintf("v_inf_earth_vector: %.4e V^ %.4e B^ km/s\n", v_inf_earth_vector);
fprintf("v_inf_earth: %.4e km/s\n", v_inf_earth);
fprintf("v_pass_earth: %.4e km/s\n", v_pass_earth);

fprintf("-------Problem 2f: Consider Mars Local Field-----------\n")
r_pass_mars = 500 + R_mars;
v_arr_vector = [ v_arr 0 ];
v_mars_vector = [ v2*cosd(-FPA_arr) v2*sind(-FPA_arr) ];
v_inf_mars_vector = v_arr_vector - v_mars_vector;
v_inf_mars = norm(v_inf_mars_vector);
energy_hyp_mars = v_inf_mars^2/2;
v_pass_mars = sqrt(2*(energy_hyp_mars + (mu_mars/r_pass_mars)));

fprintf("v_arr_vector: %.4e V^ %.4e B^ km/s\n", v_arr_vector);
fprintf("v_mars_vector: %.4e V^ %.4e B^ km/s\n", v_mars_vector);
fprintf("v_inf_mars_vector: %.4e V^ %.4e B^ km/s\n", v_inf_mars_vector);
fprintf("v_inf_mars: %.4e km/s\n", v_inf_mars);
fprintf("energy_hyp_mars: %.4e km/s\n", energy_hyp_mars);
fprintf("v_pass_mars: %.4e km/s\n", v_pass_mars);


% plot earth's orbit
p_earth = a_earth;
plot_eph(0, p_earth, 0, 0, 360)

% plot mars' orbit
p_mars = a_mars;
plot_eph(0, p_mars, 0, 0, 360)

% plot transfer
plot_eph(e, p, 0, 79.031, 139.03)


% plot function
function plot_eph(e,p,AOP, start_range, end_range)
    %orbit distance
    ta = [start_range:0.1:end_range];
    rmag = p ./ (1+e*cosd(ta-AOP));
    rplot(1,:) = rmag.*cosd(ta);
    rplot(2,:) = rmag.*sind(ta);
    
    plot(0,0,'.','MarkerSize',10,'Color','k','HandleVisibility','off'), hold on %Primary Body
    [~,indx] = min(abs(ta-AOP));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    [~,indx] = min(abs(ta-(180+AOP)));
    plot(rplot(1,indx),rplot(2,indx),'.','MarkerSize',7,'Color','k','HandleVisibility','off') %line of apsides
    plot(rplot(1,:),rplot(2,:)); %entire orbit

    title("Lambert Arc from Earth to Mars-Lillian Shido")
    xlabel("Distance [km]")
    ylabel("Distance [km]")
    
    % Keep zoom fixed
    daspect([1 1 1]); %axis equal
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]);
end