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

% Calculate the roots
y = 0; % collinear points
x=x_L1;
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

% The solution to the orbit:
t_0 = 0; % non-dim time
xi_0 = 0.01;
eta_0 = 0;
t_end = 100*pi; % non-dim time
s = (B_1 + (B_1^2 + B_2_squared)^(1/2))^(1/2);
B_3 = (s^2 + U_xx)/2*s;
orbit_results = zeros(0,3);
for t = t_0:0.1:t_end
    xi = xi_0*cos(s*(t-t_0)) + (eta_0/B_3)*sin(s*(t-t_0));
    eta = eta_0*cos(s*(t-t_0)) - B_3*xi_0*sin(s*(t-t_0));
    orbit_results(end+1,:) = [t xi eta];
end
orbit_table = array2table(orbit_results, 'VariableNames', {'time', 'xi', 'eta'});
disp(t)

% Initial velocities
xi_v_0 = eta_0*s/B_3;
eta_v_0 = -B_3*xi_0*s;
fprintf("initial velocity xi: %f\n", xi_v_0)
fprintf("initial velocity eta: %f\n", eta_v_0)

% Plot
fig1 = figure("Name", "orbit");
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 50);
orbit = plot(orbit_table, 'xi', 'eta');
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, orbit], {'Earth', 'Moon', 'L1', 'Orbit'})
title({'The Zero Velocity Curves at z=0';['in the Earth-Moon System for Jacobi Constant = ', num2str(C)]})
box on
grid on
fontsize(14, 'points')
