clear all

% Gravitational Parameters [km^3/s^2]
mu_Sun = 132712440017.99;
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
% Calculate location of Earth-Moon L1 point
% Step 1
L1_results = zeros(0,9);
L1_tolerance = 1e-12;
L1_counter = 0;
gamma = 0.2; % initial guess
while 1
    L1_counter = L1_counter + 1;
    % Step 2
    f = 1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2);
    f_prime = -1-(2*(1-mu)/(1-gamma)^3)-(2*mu/gamma^3);
    % result = 1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2);
    % Step 3
    if abs(f) > L1_tolerance
        % Step 4
        gamma = gamma - f/f_prime;
        % gamma = gamma - ((1-mu-gamma-((1-mu)/(1-gamma)^2)+(mu/gamma^2))/((-2*(1-mu)/(1-gamma)^3)-(2*mu/gamma^3)-1));
        continue
    else
        % Step 5
        x = 1 - mu - gamma;
        % Add check for accleration
        d_x = x*a+mu;
        r_x = x*a-1+mu;
        d = ((x*a+mu)^2)^(1/2);
        r = ((x*a-1+mu)^2)^(1/2);
        accel = -(1-mu)/d^3*d_x - mu/r^3*r_x;
        % Add check for partial wrt gamma
        partial = (1-mu)/(1-gamma)^2 - mu/(gamma)^2 - (1-mu-gamma);
        L1_results(end+1,:) = [x];
        L1_counter = 0;
        break
    end
end

x_L1 = L1_results(1);
% x_L1 = 8.3692e-01; % From Problem Set B
y_L1 = 0;
fprintf("L1: %f\n", x_L1)

% Location of the other Earth-Moon Libration points
x_L2 = 1.1557e+00;
y_L2 = 0;

x_L3 = -1.0051;
y_L3 = 0;

x_L4 = cosd(60) - mu;
y_L4 = sind(60);

x_L5 = cosd(60) - mu;
y_L5 = -sind(60);

% Perturbations
xi_0 = 5e-6;
xi_0_dim = xi_0*l_char;
fprintf("Dimensional xi_0: %f km\n", xi_0_dim)
eta_0 = 0;
x_0 = x_L1 + xi_0;
y_0 = 0 + eta_0;

% The solution to the perturbed orbit:
t_0 = 0; % non-dim time
t_end = 1*pi; % non-dim time
orbit_results = zeros(0,4);
% d, r, and partials all evaluated @ L1
d = sqrt((x_L1+mu)^2 + y_L1^2);
r = sqrt((x_L1-1+mu)^2 + y_L1^2);
U_xx = 1 - (1-mu)/d^3 - mu/r^3 + 3*(1-mu)*(x_L1+mu)^2/d^5 + 3*mu*(x_L1-1+mu)^2/r^5;
fprintf("Uxx: %f\n", U_xx)
U_yy = 1 - (1-mu)/d^3 - mu/r^3;
B_1 = 2 - (U_xx + U_yy)/2;
B_2_squared = -U_xx*U_yy;
s = (B_1 + (B_1^2 + B_2_squared)^(1/2))^(1/2);
B_3 = (s^2 + U_xx)/(2*s);
big_Lambda_1 = -B_1 + (B_1^2 + B_2_squared)^(1/2);
big_Lambda_2 = -B_1 - (B_1^2 + B_2_squared)^(1/2);
lambda_1 = sqrt(big_Lambda_1); % real
lambda_2 = -sqrt(big_Lambda_1); % real
lambda_3 = sqrt(big_Lambda_2); % imaginary
lambda_4 = -sqrt(big_Lambda_2); % imaginary
for t = t_0:0.1:t_end
    xi = xi_0*cos(s*(t-t_0)) + (eta_0/B_3)*sin(s*(t-t_0));
    x = xi + x_L1;
    eta = eta_0*cos(s*(t-t_0)) - B_3*xi_0*sin(s*(t-t_0)); 
    y = eta + y_L1; 
    orbit_results(end+1,:) = [t x y x-x_L1];
end
orbit_table = array2table(orbit_results, 'VariableNames', {'time', 'x', 'y', 'xi'});
disp(t)

% Calculate the roots of the L1 point
fprintf("x_L: %f\n", x)
fprintf("lambda_1: %f + %fi\n", real(lambda_1), imag(lambda_1))
fprintf("lambda_2: %f + %fi\n", real(lambda_2), imag(lambda_2))
fprintf("lambda_3: %f + %fi\n", real(lambda_3), imag(lambda_3))
fprintf("lambda_4: %f + %fi\n", real(lambda_4), imag(lambda_4))

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

% Calculate the Jacobi Constant for these ICs. 
x_0 = x_L1 + xi_0;
y_0 = 0;
d_0 = sqrt((x_0+mu)^2 + y_0^2);
r_0 = sqrt((x_0-1+mu)^2 + y_0^2);
pseudo_U = (1-mu)/d_0 + mu/r_0 + (x_0^2+y_0^2)/2;
v_squared = v_x0^2 + v_y0^2;
C = 2*pseudo_U - v_squared;
fprintf("d_0: %f\n", d_0)
fprintf("r_0: %f\n", r_0)
fprintf("pseudo_U: %f\n", pseudo_U)
fprintf("v_squared: %f\n", v_squared)
fprintf("Jacobi Constant C: %d\n", C)

% To plot the ZVC for z = 0:
% Step 1: For a given x
% Step 2: Guess a y
% Step 3: Calculate Jacobi constant result
% Step 4: If Jacobi constant result is less than tolerance, break
% Step 5: If not, update y with Newton Raphson
% Step 6: Repeat

zvc_result = zeros(0,7);
tolerance = 1e-12;
max_iterations = 1000;

%===================Calculate the ZVC around the Earth=====================
%=============Useful Functions=========================

% Find Jacobi Constant C
function C = find_JC(x, y, mu)
    d = sqrt((x+mu)^2 + y^2);
    r = sqrt((x-1+mu)^2 + y^2);
    C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
end

% Calc error
function error = calc_error(actual, ideal)
    error =  abs(actual - ideal)/ideal;
end

% Give an initial guess for y that is around the Earth)
% for x = linspace(-1.5,1.5,5e3) %  Find the curve for -1.5 < x < 1.5
%     for y = 0:0.1:1.5 % These are my guesses
%         counter = 0;
%         while 1
%             counter = counter + 1;
%             d = sqrt((x+mu)^2 + y^2);
%             r = sqrt((x-1+mu)^2 + y^2);
%             f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%             f_prime_y = 2*y*( 1 - (1-mu)/d^3 - mu/r^3);
%             delta = f/f_prime_y;
%             if abs(f) > tolerance
%                 if counter > max_iterations % check for max iterations
%                     new_C = find_JC(x, y, mu);
%                     C_error = calc_error(new_C, C);
%                     if C_error < 1e-12
%                         zvc_result(end+1,:) = [x y -y x-x_L1 new_C error_C counter];
%                     else
%                         break
%                     end
%                 elseif d == 0 || r == 0 % check for dividing by 0
%                     break
%                 elseif abs(f_prime_y) <= tolerance % check for zero derivative
%                     break
%                 else
%                     y = y - delta;
%                     continue
%                 end
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y -y x-x_L1 new_C error_C counter];
%                 break
%             end
%         end
%     end
% end
% for y = linspace(-1.5,1.5,5e3) % Find the curve for 0 < y < 1.21
%     for x = -1.5:0.1:1.5 % These are my guesses for x
%         counter = 0;
%         while 1
%            counter = counter + 1;
%            d = sqrt((x+mu)^2 + y^2);
%            r = sqrt((x-1+mu)^2 + y^2);
%            f = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r) - C;
%            f_prime_x = 2*x - 2*(1-mu)*(x+mu)/d^3 - 2*mu*(x-1+mu)/r^3;
%            delta = f/f_prime_x;
%            if abs(f) > tolerance
%                 if counter > max_iterations % check for max iterations
%                     new_C = find_JC(x, y, mu);
%                     C_error = calc_error(new_C, C);
%                     if C_error < 1e-12
%                         zvc_result(end+1,:) = [x y -y x-x_L1 new_C error_C counter];
%                     else
%                         break
%                     end
%                 elseif d == 0 || r == 0 % check for dividing by 0
%                     break
%                 elseif abs(f_prime_x) <= tolerance % check for zero derivative
%                     break
%                 else
%                     x = x - delta;
%                     continue
%                 end
%             else
%                 new_C = x^2 + y^2 + (2*(1-mu)/d) + (2*mu/r);
%                 error_C = abs(new_C - C)/C*100;
%                 zvc_result(end+1,:) = [x y -y x-x_L1 new_C error_C counter];
%                 break
%             end
%         end
%     end
% end

zvc_table = array2table(zvc_result, 'VariableNames', {'x','y', '-y', 'xi', 'JC', '% Error of C','Iterations'});
format short
disp(zvc_table)

% Find the max error of JC
max_error = max(abs(zvc_result(:,4)));
fprintf("max error: %d\n", max_error)

%=================Calculate the non-linear result=======================

% Initial Conditions
w0 = [x_0;y_0;v_x0;v_y0];

% tspan should be calculated with non-dim time
tspan = [0 t_end];

% Setting up the ODE
ode = @(t,w) [...
w(3);...
w(4);...
2*w(4) + w(1) - (1 - mu) * (w(1) + mu) / ((w(1) + mu)^2 + w(2)^2)^(3/2)...
- mu * (w(1) - 1 + mu) / ((w(1) - 1 + mu)^2 + w(2)^2)^(3/2);...
-2*w(3) + w(2) - (1 - mu) * w(2) / ((w(1) + mu)^2 + w(2)^2)^(3/2) - mu * w(2)/((w(1) - 1 + mu)^2 + w(2)^2)^(3/2)...
];
options = odeset('RelTol',1e-12,'AbsTol', 1e-14,'MaxStep', 1e-3);
[t,w] = ode113(ode, tspan, w0, options); 

% Retrieving the results from the integrator
% Position (dimensional)
x_nl = w(:,1);
y_nl = w(:,2);
% Position (non-dimensional)
xi_nl = w(:,1)-x_L1;
eta_nl = w(:,2)-0;
% Velocity
v_x_nl = w(:,3);
v_y_nl = w(:,4);

% Calculate the Jacobi Constant C
d = sqrt((x_nl+mu).^2 + y_nl.^2);
r = sqrt((x_nl-1+mu).^2 + y_nl.^2);
x_y_sq = (x_nl.^2+y_nl.^2)./2;
term_1 = (1-mu)./d;
term_2 = mu./r;
pseudo_U = term_1 + term_2 + x_y_sq;
v_squared = w(:,3).^2 + w(:,4).^2;
C_table = 2*pseudo_U - v_squared;
error_table = abs(C_table - C)./C; 

%===================Plot the orbit===============
fig1 = figure("Name", "orbit");
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 10);
linear_orbit = plot(orbit_table, 'x', 'y');
nonlinear_orbit = plot(x_nl, y_nl);
hold off
limit = 5*xi_0;
xlim([-limit limit])
ylim([-limit limit])
axis square
xlabel("Non-dimensional X")
ylabel("Non-dimensional Y")
legend([earth, moon, L1_plot, linear_orbit, nonlinear_orbit], {'Earth', 'Moon', 'L1', 'Linear Orbit', 'Non-Linear Orbit'})
title({'The orbit around the L1 point in the Earth-Moon System'})
box on
grid on
fontsize(14, 'points')

%===================Plot the ZVC!===============%
fig1 = figure('Name','ZVC');
earth = scatter(x_Earth, 0, 'blue', 'filled', 'SizeData', 200);
hold on
moon = scatter(x_Moon, 0, 'black', 'filled', 'SizeData', 100);
L1_plot = scatter(x_L1, y_L1, 'red', 'filled', 'SizeData', 20);
scatter(x_L2, y_L2, 'red', 'filled', 'SizeData', 20)
scatter(x_L3, y_L3, 'red', 'filled', 'SizeData', 20)
scatter(x_L4, y_L4, 'red', 'filled', 'SizeData', 20)
scatter(x_L5, y_L5, 'red', 'filled', 'SizeData', 20)
zvc_plot = scatter(zvc_table, 'x', 'y', 'filled', 'SizeData', 5, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9');
scatter(zvc_table, 'x', '-y', 'filled', 'SizeData', 5, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9')
orbit = plot(orbit_table, 'x', 'y');
hold off
xlim([-1.5 1.5])
ylim([-1.5 1.5])
axis square
xlabel("X [non-dim]")
ylabel("Y [non-dim]")
legend([earth, moon, L1_plot, zvc_plot, orbit], {'Earth', 'Moon', 'Eq. Points', 'ZVC', 'Orbit'})
title({'The Zero Velocity Curves at z=0';['in the Earth-Moon System for Jacobi Constant = ', num2str(C)]})
box on
grid on
fontsize(14, 'points')

fig2 = figure('Name','Error C');
error_C = scatter(t, error_table,'SizeData', 5);
hold on
hold off
xlabel("Time [rad]")
ylabel("Jacobi Constant Error")
title('Error of the Jacobi Constant over time')
box on
grid on
fontsize(14, 'points')

fig3 = figure('Name', 'non-dim');
L1_plot = scatter(0, 0, 'red', 'filled', 'SizeData', 10);
hold on
% zvc_plot = scatter(zvc_table, 'xi', 'y', 'filled', 'SizeData', 2, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9');
% scatter(zvc_table, 'xi', '-y', 'filled', 'SizeData', 2, 'MarkerFaceColor','#53A1C9', 'MarkerEdgeColor', '#53A1C9')
linear_orbit = plot(orbit_table, 'xi', 'y', 'Color', '#800080');
nonlinear_orbit = plot(xi_nl, eta_nl, 'Color', '#008000');
hold off
xlim([-limit limit])
ylim([-limit limit])
axis square
xlabel("\xi [non-dim]")
ylabel("\eta [non-dim]")
legend([L1_plot, zvc_plot, linear_orbit, nonlinear_orbit], {'L1', 'ZVC', 'Linear Orbit', 'Non-Linear Orbit'})
title({'The orbit around the L1 point in'; ['the Earth-Moon System for \xi = ', num2str(xi_0)]})
box on
grid on
fontsize(14, 'points')