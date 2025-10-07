clear all

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
mu_Earth = 398600.4415;
mu_Moon = 4902.8005821478;

% Earth-Moon System
mu = mu_Moon/(mu_Earth + mu_Moon);
% Position of primary bodies
x_Earth = -mu;
x_Moon = 1-mu;
% Location of Earth-Moon L1 point
x_L1 = 8.3692e-01; % From Problem Set B
y_L1 = 0;
% Characteristic Length [km]
a_Moon = 384400; % around Earth
l_char = a_Moon;
fprintf("Characteristic Length: %f km \n", l_char)
% Calculate characteristic time
t_char = sqrt(a_Moon^3/(mu_Earth+mu_Moon));
fprintf("characteristic time: %d sec\n", t_char)

% Calculate the roots of the L1 point
y = 0; % collinear points
x = x_L1;
d = sqrt((x+mu)^2 + y^2);
r = sqrt((x-1+mu)^2 + y^2);
U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x+mu)^2/d^5 + 3*mu*(x-1+mu)^2/r^5;
U_yy = 1 - (1-mu)/d^3 - mu/r^3;
B_1 = 2 - (U_xx + U_yy)/2;
B_2_squared = -U_xx*U_yy;
big_Lambda_1 = -B_1 + (B_1^2 + B_2_squared)^(1/2);
big_Lambda_2 = -B_1 - (B_1^2 + B_2_squared)^(1/2);
lambda_1 = sqrt(big_Lambda_1); % real
lambda_2 = -sqrt(big_Lambda_1); % real
lambda_3 = sqrt(big_Lambda_2); % imaginary
lambda_4 = -sqrt(big_Lambda_2); % imaginary
fprintf("x_L: %f\n", x)
fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))

% Perturbations
xi_0 = 0.01;
xi_0_dim = xi_0*l_char;
fprintf("Dimensional xi_0: %f km\n", xi_0_dim)
eta_0 = 0;
x_0 = x_L1 + xi_0;
y_0 = 0 + eta_0;

% The solution to the perturbed orbit:
t_0 = 0; % non-dim time
t_end = 8*pi; % non-dim time
orbit_results = zeros(0,5);
% d, r, and partials all evaluated @ L1
d = sqrt((x_L1+mu)^2 + y_L1^2);
r = sqrt((x_L1-1+mu)^2 + y_L1^2);
U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x_L1+mu)^2/d^5 + 3*mu*(x_L1-1+mu)^2/r^5;
U_yy = 1 - (1-mu)/d^3 - mu/r^3;
B_1 = 2 - (U_xx + U_yy)/2;
B_2_squared = -U_xx*U_yy;
s = (B_1 + (B_1^2 + B_2_squared)^(1/2))^(1/2);
B_3 = (s^2 + U_xx)/(2*s);
for t = t_0:0.1:t_end
    xi = xi_0*cos(s*(t-t_0)) + (eta_0/B_3)*sin(s*(t-t_0));
    x = xi + x_L1;
    eta = eta_0*cos(s*(t-t_0)) - B_3*xi_0*sin(s*(t-t_0)); 
    y = eta + y_L1; 
    orbit_results(end+1,:) = [t x y xi eta];
end
orbit_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y', 'xi', 'eta'});
disp(t)

% Max distance from L1
% Where y is maximum
y_max = eta_0*cosd(90) - B_3*xi_0*sind(90) + y_L1;
fprintf('Max distance from L1 (non-dim): %f\n', y_max)
fprintf('Max distance from L1 (dim): %f\n', y_max*l_char)

% Period of the orbit around L1
P_around_L1 = 2*pi/s;
fprintf('Period around L1 (non-dim relative to the larger system): %f\n', P_around_L1)
fprintf('Period around L1 (dim): %f sec\n', P_around_L1*t_char)
fprintf('Period around L1 (dim): %f days\n', P_around_L1*t_char/3600/24)

% Initial velocities
v_x0 = eta_0*s/B_3;
v_y0 = -B_3*xi_0*s;
fprintf("initial velocity x-comp: %f\n", v_x0)
fprintf("initial velocity y-comp: %f\n", v_y0)
v_x0_dim = v_x0 * l_char / t_char * 1000;
v_y0_dim = v_y0 * l_char / t_char * 1000;
fprintf("initial velocity x-comp (dim): %f m/s\n", v_x0_dim)
fprintf("initial velocity y-comp (dim): %f m/s\n", v_y0_dim)

% Plot
fig1 = figure("Name", "orbit");
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 10);
orbit = plot(orbit_table, 'x', 'y');
hold off
xlim([-0.2 1.2])
ylim([-0.7 0.7])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, orbit], {'Earth', 'Moon', 'L1', 'Orbit'})
title({'The orbit around the L1 point in the Earth-Moon System'})
box on
grid on
fontsize(14, 'points')

fig2 = figure("Name", "xi and eta")
L1_plot = scatter(0, 0, 'red', 'filled', 'SizeData', 20);
hold on
orbit = plot(orbit_table, 'xi', 'eta');
xlim([-0.05 0.05])
ylim([-0.05 0.05])
axis square
grid on
box on
title({'The orbit relative to the L1 point in the Earth-Moon System'})
xlabel("Non-dimensional \xi")
ylabel("Non-dimensional \eta")
legend([L1_plot, orbit], {'L1', 'Orbit'})
fontsize(14, 'points')
