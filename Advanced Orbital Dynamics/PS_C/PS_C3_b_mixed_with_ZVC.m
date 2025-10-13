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

% Calculate the linear orbit for the equilateral points:
function linear_table = calc_equilateral_linear(t0, t_end, alpha_1, alpha_2, alpha_3, alpha_4, beta_1, beta_2, beta_3, beta_4, s1, s2, x_L, y_L)
    orbit_results = zeros(0,5);
    for t = t0:0.01:t_end
        xi = alpha_1*cos(s1*t) + alpha_2*sin(s1*t) + alpha_3*cos(s2*t) + alpha_4*sin(s2*t);
        eta = beta_1*cos(s1*t) + beta_2*sin(s1*t) + beta_3*cos(s2*t) + beta_4*sin(s2*t);
        x = xi + x_L;
        y = eta + y_L; 
        orbit_results(end+1,:) = [t x y xi eta];
    end
    linear_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'});
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

% Calculate alpha_2 and beta_2 for equilateral points:
function [alpha_2, beta_2] = calc_alpha2_beta2(U_xx, U_yy, U_xy, xi_0, eta_0, s1)
    alpha_2 = (U_xy*xi_0 + U_yy*eta_0 + eta_0*s1^2)/(2*s1);
    beta_2 = (U_xx*xi_0 + U_xy*eta_0 + xi_0*s1^2)/(-2*s1);
end

% Calculate alpha_4 and beta_4 for equilateral points:
function [alpha_4, beta_4] = calc_alpha4_beta4(U_xx, U_yy, U_xy, xi_0, eta_0, s2)
    alpha_4 = (U_xy*xi_0 + U_yy*eta_0 + eta_0*s2^2)/(2*s2);
    beta_4 = (U_xx*xi_0 + U_xy*eta_0 + xi_0*s2^2)/(-2*s2);
end

% Calculate the Jacobi Constant
function C = calc_JC(mu, x, y, x_dot, y_dot)
    d = sqrt((x+mu)^2 + y^2);
    r = sqrt((x-1+mu)^2 + y^2);
    pseudo_U = (1-mu)/d + mu/r + (x^2+y^2)/2;
    v_squared = (x_dot^2 + y_dot^2);
    C = 2*pseudo_U - v_squared;
end

% Find Jacobi Constant C
function C = JC_check(mu, x, y)
    d = sqrt((x+mu)^2 + y^2);
    r = sqrt((x-1+mu)^2 + y^2);
    C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
end

% Calc error
function error = calc_error(actual, ideal)
    error =  abs(actual - ideal)/ideal;
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
% Find a set of ICs that include both short and long periods for the motion. 
% Maybe we should take the long and short period ICs and add them?
% Perturbations
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

% Calc s
s1 = imag(lambda_1);
s2 = imag(lambda_3);
fprintf("s1: %f, s2: %f\n", s1, s2)
fprintf("s1/s2 ratio: %f\n", s1/s2)
fprintf("s2/s1 ratio: %f\n", s2/s1)

% Set combined alphas:
% alpha_1 = xi_0;
% beta_1 = eta_0;
% [alpha_2, beta_2] = calc_alpha2_beta2(U_xx, U_yy, U_xy, xi_0, eta_0, s1);
% alpha_3 = xi_0;
% beta_3 = eta_0;
% [alpha_4, beta_4] = calc_alpha4_beta4(U_xx, U_yy, U_xy, xi_0, eta_0, s2);
% Long Period
alpha_1 = 0.1;
beta_1 = 0;
alpha_2 = 0.21251;
beta_2 = -0.14066;
% Short Period
alpha_3 = 0.01;
beta_3 = 0;
alpha_4 = 0.006639;
beta_4 = -0.008701;
fprintf("alpha_1: %f, beta_1: %f\n", alpha_1, beta_1)
fprintf("alpha_2: %f, beta_2: %f\n", alpha_2, beta_2)
fprintf("alpha_3: %f, beta_3: %f\n", alpha_3, beta_3)
fprintf("alpha_4: %f, beta_4: %f\n", alpha_4, beta_4)

t0 = 0;
t_end = 50*pi;
linear_table = calc_equilateral_linear(t0, t_end, alpha_1, alpha_2, alpha_3, alpha_4, beta_1, beta_2, beta_3, beta_4, s1, s2, x_L4, y_L4);
% orbit_table: 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'}

% Calc initial conditions now
xi_0 = alpha_1 + alpha_3;
eta_0 = beta_1 + beta_3;
xi_dot_0 = alpha_2*s1 + alpha_4*s2;
eta_dot_0 = beta_2*s1 + beta_4*s2;
fprintf("xi_0: %f, eta_0: %f [non-dim]\n", xi_0, eta_0)
fprintf("xi_0: %f, eta_0: %f [km]\n", xi_0*l_char, eta_0*l_char)
fprintf("xi_dot_0: %f, eta_dot_0: %f [non-dim]\n", xi_dot_0, eta_dot_0)
fprintf("xi_dot_0: %f, eta_dot_0: %f [m/s]\n", xi_dot_0*l_char/t_char*1000, eta_dot_0*l_char/t_char*1000)

% Calc Jacobi Constant
C = calc_JC(mu, x_L4+xi_0, y_L4+eta_0, xi_dot_0, eta_dot_0);
% C = 3;
fprintf("Jacobi Constant: %f\n", C)

%===========================END DEFINE PROBLEM=============================

%==============================CALC ZVC===================================
% Calc ZVCs
zvc_result = zeros(0,6);
tolerance = 1e-12;
max_iterations = 100;
% Give an initial guess for y that is around the Earth
% for x = linspace(-1.5,1.5,1e3) %  Find the curve for -1.5 < x < 1.5
for x = linspace(0,1,5e2) %  Find the curve for -1.5 < x < 1.5
    fprintf("x sweep: %f\n", x)
    for y = 0:0.1:1.5 % These are my guesses
        counter = 0;
        while 1
            counter = counter + 1;
            d = sqrt((x+mu)^2 + y^2);
            r = sqrt((x-1+mu)^2 + y^2);
            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
            f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
            delta = f/f_prime_y;
            if abs(f) > tolerance
                if counter > max_iterations % check for max iterations
                    new_C = JC_check(mu, x, y);
                    C_error = calc_error(new_C, C);
                    if C_error < 1e-12
                        zvc_result(end+1,:) = [x y -y new_C error_C counter];
                    else
                        break
                    end
                elseif d == 0 || r == 0 % check for dividing by 0
                    break
                elseif abs(f_prime_y) <= tolerance % check for zero derivative
                    break
                else
                    y = y - delta;
                    continue
                end
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                break
            end
        end
    end
end
% for y = linspace(-1.5,1.5,1e3) % Find the curve for 0 < y < 1.21
for y = linspace(0.5,1.1,5e2) % Find the curve for 0 < y < 1.21
    fprintf("y sweep: %f\n", y)
    if y == 0.752505
        fprintf("y = 0.752505!!\n")
    end
    for x = -1.5:0.1:1.5 % These are my guesses for x
        counter = 0;
        while 1
        counter = counter + 1;
        d = sqrt((x+mu)^2 + y^2);
        r = sqrt((x-1+mu)^2 + y^2);
        f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
        f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
        delta = f/f_prime_x;
        if abs(f) > tolerance
                if counter > max_iterations % check for max iterations
                    new_C = JC_check(mu, x, y);
                    C_error = calc_error(new_C, C);
                    if C_error < 1e-12
                        zvc_result(end+1,:) = [x y -y new_C error_C counter];
                    else
                        break
                    end
                elseif d == 0 || r == 0 % check for dividing by 0
                    break
                elseif abs(f_prime_x) <= tolerance % check for zero derivative
                    break
                else
                    x = x - delta;
                    continue
                end
            else
                new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
                error_C = abs(new_C - C)/C*100;
                zvc_result(end+1,:) = [x y -y new_C error_C counter];
                break
            end
        end
    end
end
zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', '-y', 'JC', '% Error of C','Iterations'});
%============================END CALC ZVC==================================

%===========================CONFIGURE PLOTS================================
fig1 = figure('Name', '1 :non-dim');
L4_plot = scatter(0, 0, 'red', 'filled', 'SizeData', 10);
hold on
linear_orbit = plot(linear_table, 'xi', 'eta', 'Color', '#800080');
hold off
% limit = 5*xi_0;
% ylim([-limit limit])
% xlim([-limit limit])
axis square
xlabel("\xi [non-dim]")
ylabel("\eta [non-dim]")
% xticks(-limit:.01:limit)
legend([L4_plot, linear_orbit], {'L4', 'Linear Orbit'})
title({'L4 Mixed frequency linear result (\xi-\eta) in'; 'the Earth-Moon System, t=50\pi (Lillian Shido)'})
box on
grid on
fontsize(14, 'points')

fig2 = figure('Name', '2: dim');
L4_plot = scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 10);
hold on
% earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
% moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
center = scatter(0, 0);
linear_orbit = plot(linear_table, 'x', 'y', 'Color', '#800080');
zvc = scatter(zvc_table, 'x', 'y', 'SizeData', 2, 'MarkerFaceColor', '#53A1C9', 'MarkerEdgeColor', '#53A1C9');
% scatter(zvc_table, 'x', '-y', 'SizeData', 5, 'MarkerFaceColor', '#53A1C9', 'MarkerEdgeColor', '#53A1C9')
hold off
xlim([0 0.8])
ylim([0.5 1.3])
axis square
xlabel("x [non-dim]")
ylabel("y [non-dim]")
legend([L4_plot, linear_orbit, zvc], {'L4', 'Linear Orbit', 'ZVC'})
title({'L4 Mixed frequency linear (x-y)'; ['JC=', num2str(C),', t=',num2str(t_end),'rad (Lillian Shido)']})
box on
grid on
fontsize(14, 'points')