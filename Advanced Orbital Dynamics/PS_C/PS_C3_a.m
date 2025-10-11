clear all

%=========================FUNCTION DEFINITIONS============================%

% Calculate the location of L4
function [x_L4, y_L4] = calc_L4(mu)
    x_L4 = cosd(60) - mu;
    y_L4 = sind(60);
end

% Calculate Uxx, Uyy, Uxy(=Uyx)
function [U_xx, U_yy, U_xy] = calc_U(mu, x_L, y_L)
    d = sqrt((x_L+mu)^2 + y_L^2);
    r = sqrt((x_L-1+mu)^2 + y_L^2);
    U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x_L+mu)^2/d^5 + 3*mu*(x_L-1+mu)^2/r^5;
    U_yy = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*y_L^2/d^5 + 3*mu*y_L^2/r^5;
    U_xy = 3*(1-mu)*(x_L+mu)*y_L/d^5 + 3*mu*(x_L-1+mu)*y_L/r^5;
end

% Calculate eigenvalues for equilateral points:
function [lambda_1, lambda_2, lambda_3, lambda_4] = calc_eigenvalues_equilateral(mu)
    big_lambda_1 = 1/2*(-1 + (1 - 27*mu*(1-mu))^(1/2));
    big_lambda_2 = 1/2*(-1 - (1 - 27*mu*(1-mu))^(1/2));
    lambda_1 = sqrt(big_lambda_1);
    lambda_2 = -lambda_1;
    lambda_3 = sqrt(big_lambda_2);
    lambda_4 = -lambda_3;
end

% Calculate alphas for equilateral points:
function [alpha, check_alpha] = calc_alpha(lambda, U_xx, U_yy, U_xy)
    alpha = (lambda^2-U_xx)/(2*lambda + U_xy);
    check_alpha = -(2*lambda - U_xy)/(lambda^2 - U_yy);
end

% Calculate the long period linear orbit for the equilateral points:
function linear_table = calc_equilateral_linear_long(t0, t_end, alpha_1, alpha_2, beta_1, beta_2, s1, x_L, y_L)
    orbit_results = zeros(0,5);
    for t = t0:0.01:t_end
        xi = alpha_1*cos(s1*t) + alpha_2*sin(s1*t);
        eta = beta_1*cos(s1*t) + beta_2*sin(s1*t);
        x = xi + x_L;
        y = eta + y_L; 
        orbit_results(end+1,:) = [t x y xi eta];
    end
    linear_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'});
end

%===========================END FUNCTION DEFINITIONS===============================

%==========================DEFINE SYSTEM==================================

% Earth-Moon System
% Gravitational Parameters [km^3/s^2]
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Characteristic Length [km]
a = 384400; % around Earth
l_char = a;
fprintf("Characteristic Length: %f km \n", l_char)
% Calculate characteristic time
t_char = sqrt(a^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% Earth-Moon System
mu = mu_Moon/(mu_Earth + mu_Moon);
fprintf("Mu of Earth-Moon: %f\n", mu)
% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;

% Calc L4 Libration point
[x_L4, y_L4] = calc_L4(mu);
fprintf("x_L4: %f, y_L4: %f\n", x_L4, y_L4)

%============================END DEFINE SYSTEM================================
%===========================DEFINE PROBLEM=============================
% Define Perturbations
xi_0 = 0.01;
eta_0 = 0;

% Calc Partials of Pseudo-potential
[U_xx, U_yy, U_xy] = calc_U(mu, x_L4, y_L4);
fprintf("U_xx: %f\n", U_xx)
fprintf("U_yy: %f\n", U_yy)
fprintf("U_xy: %f\n", U_xy)

% Calc eigenvalues for L4
[lambda_1, lambda_2, lambda_3, lambda_4] = calc_eigenvalues_equilateral(mu);
fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))

% Calc alphas for L4
[alpha_1, check_alpha_1] = calc_alpha(lambda_1, U_xx, U_yy, U_xy);
fprintf("alpha_1: %f + %fi\n", real(alpha_1), imag(alpha_1))
[alpha_2, check_alpha_2] = calc_alpha(lambda_2, U_xx, U_yy, U_xy);
fprintf("alpha_2: %f + %fi\n", real(alpha_2), imag(alpha_2))
[alpha_3, check_alpha_3] = calc_alpha(lambda_3, U_xx, U_yy, U_xy);
fprintf("alpha_3: %f + %fi\n", real(alpha_3), imag(alpha_3))
[alpha_4, check_alpha_4] = calc_alpha(lambda_4, U_xx, U_yy, U_xy);
fprintf("alpha_4: %f + %fi\n", real(alpha_4), imag(alpha_4))

% Calc betas for L4
% beta_1 = alpha_1/i;
% beta_2 = alpha_2/i;
% beta_3 = alpha_3/i;
% beta_4 = alpha_4/i;
beta_1 = imag(alpha_1);
beta_2 = imag(alpha_2);
beta_3 = imag(alpha_3);
beta_4 = imag(alpha_4);
fprintf("beta_1: %f + %fi\n", real(beta_1), imag(beta_1))
fprintf("beta_2: %f + %fi\n", real(beta_2), imag(beta_2))
fprintf("beta_3: %f + %fi\n", real(beta_3), imag(beta_3))
fprintf("beta_4: %f + %fi\n", real(beta_4), imag(beta_4))

% Calc s
s1 = imag(lambda_1);
s2 = imag(lambda_3);
fprintf("s1: %f, s2: %f\n", s1, s2)

% Calc the longer period linear orbit
t0 = 0;
t_end = 4*pi;
linear_table = calc_equilateral_linear_long(t0, t_end, alpha_1, alpha_2, beta_1, beta_2, s1, x_L4, y_L4);
% orbit_table: 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'}
%===========================END DEFINE PROBLEM=============================

%===========================CONFIGURE PLOTS================================
fig4 = figure('Name', '1 :non-dim');
L4_plot = scatter(0, 0, 'red', 'filled', 'SizeData', 10);
hold on
linear_orbit = plot(linear_table, 'xi', 'eta', 'Color', '#800080');
hold off
limit = 5*xi_0;
ylim([-limit limit])
xlim([-limit limit])
axis square
xlabel("\xi [non-dim]")
ylabel("\eta [non-dim]")
legend([L4_plot, moon, linear_orbit], {'L4', 'Moon', 'Linear Orbit'})
title({'The orbit around the L4 point in'; ['the Earth-Moon System for \xi = ', num2str(xi_0)]})
box on
grid on
fontsize(14, 'points')